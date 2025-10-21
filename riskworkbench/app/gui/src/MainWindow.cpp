#include "MainWindow.hpp"
#include "ui_MainWindow.h"

#include <QMessageBox>
#include <QTableWidgetItem>
#include <QMetaObject>
#include <QMetaType>
#include <QStatusBar>
#include <QSignalBlocker>
#include <QDebug>


#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QBarCategoryAxis>

#include <QHBoxLayout>
#include <QHeaderView>
#include <QTableWidget>
#include <QProgressBar>
#include <QComboBox>
#include <QLabel>
#include <QPushButton>
#include <QImage>
#include <QPixmap>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QPainter>

#include <QFileDialog>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDir>
#include <QDateTime>
#include <QFileInfo>

#include <QShortcut>
#include <QKeySequence>
#include <QTabWidget>

#include <QDesktopServices>
#include <QUrl>
#include <QClipboard>
#include <QGuiApplication>

#include <cmath>               
#include <limits>  
#include <algorithm>
#include <numeric>

#include "McWorker.hpp"
#include "CalibWorker.hpp"

// RW
#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/qc/qc.hpp>

namespace {
const QColor C_MEAN (255,193,7);     // MC mean (amber)
const QColor C_CI   (33,150,243);    // CI band (blue)
const int    A_CI   = 64;            // alpha for CI fill
const QColor C_BS   (76,175,80);     // BS ref = GREEN SOLID
const QColor C_LAST (156,39,176);    // Last batch point (purple)
}

namespace {

// CDF/PDF normales
inline double norm_cdf(double x) { return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0))); }
inline double norm_pdf(double x) { static const double INV_SQRT2PI = 0.3989422804014327;
                                   return INV_SQRT2PI * std::exp(-0.5*x*x); }

// ---- PRIX BS (retourne un double) ----
inline double bs_price(const rw::market::MarketData& mkt,
                       const rw::market::Instrument& inst,
                       double sigma)
{
  const double S = mkt.S0, K = inst.K, r = mkt.r, q = mkt.q, T = std::max(0.0, inst.T);
  if (T == 0.0) {
    return (inst.type == rw::market::OptionType::Call) ? std::max(S - K, 0.0)
                                                       : std::max(K - S, 0.0);
  }
  const double vol = std::max(1e-16, sigma);
  const double sqrtT = std::sqrt(T);
  const double d1  = (std::log(S/K) + (r - q + 0.5*vol*vol)*T) / (vol*sqrtT);
  const double d2  = d1 - vol*sqrtT;
  const double dfq = std::exp(-q*T);
  const double dfr = std::exp(-r*T);

  if (inst.type == rw::market::OptionType::Call)
    return S*dfq*norm_cdf(d1) - K*dfr*norm_cdf(d2);
  else
    return K*dfr*norm_cdf(-d2) - S*dfq*norm_cdf(-d1);
}

// ---- GRECS BS (pour le test Greeks vs BS) ----
struct BsGreeks { double delta{}, gamma{}, vega{}, rho{}, theta{}; };

inline BsGreeks bs_greeks(const rw::market::MarketData& mkt,
                          const rw::market::Instrument& inst,
                          double sigma)
{
  BsGreeks g{};
  const double S = mkt.S0, K = inst.K, r = mkt.r, q = mkt.q, T = std::max(0.0, inst.T);
  if (T == 0.0) {
    const bool call = (inst.type == rw::market::OptionType::Call);
    g.delta = call ? (S > K ? 1.0 : 0.0) : (S < K ? -1.0 : 0.0);
    // gamma, vega, rho, theta ~ 0 à T=0 (hors frontière)
    return g;
  }

  const double vol   = std::max(1e-16, sigma);
  const double sqrtT = std::sqrt(T);
  const double d1    = (std::log(S/K) + (r - q + 0.5*vol*vol)*T) / (vol*sqrtT);
  const double d2    = d1 - vol*sqrtT;
  const bool call    = (inst.type == rw::market::OptionType::Call);
  const double dfq   = std::exp(-q*T);
  const double dfr   = std::exp(-r*T);

  const double Nd1   = norm_cdf(d1), Nd2 = norm_cdf(d2);
  const double Nmd1  = norm_cdf(-d1), Nmd2 = norm_cdf(-d2);
  const double phid1 = norm_pdf(d1);

  g.delta = call ? dfq*Nd1 : dfq*(Nd1 - 1.0);
  g.gamma = dfq*phid1 / (S*vol*sqrtT);
  g.vega  = S*dfq*phid1*sqrtT;                 // par 1 pt de sigma (absolu)
  g.rho   = call ?  K*T*dfr*Nd2 : -K*T*dfr*Nmd2;
  const double term = S*dfq*phid1*vol / (2.0*sqrtT);
  g.theta = call ? (-term - r*K*dfr*Nd2 + q*S*dfq*Nd1)
                 : (-term + r*K*dfr*Nmd2 - q*S*dfq*Nmd1);
  return g;
}

} // namespace

static QJsonValue labelToJson(const QLabel* L) {
  if (!L) return QJsonValue(); // null
  bool ok=false; const double v = L->text().toDouble(&ok);
  return ok ? QJsonValue(v) : QJsonValue(); // null si vide/non-numérique
}

static void setLabelFromJson(QLabel* L, const QJsonObject& o, const char* key, int prec=6) {
  if (!L) return;
  if (!o.contains(key)) return;
  const QJsonValue v = o.value(key);
  if (v.isDouble()) L->setText(QString::asprintf(("%." + QString::number(prec) + "f").toStdString().c_str(),
                                                 v.toDouble()));
  else if (v.isString()) L->setText(v.toString());
}

MainWindow::MainWindow(QWidget* parent)
  : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  // Enregistrements pour les queued connections
  qRegisterMetaType<std::size_t>("std::size_t");
  qRegisterMetaType<rw::calib::CalibReport>("rw::calib::CalibReport");
  qRegisterMetaType<std::vector<rw::qc::ErrorRow>>("std::vector<rw::qc::ErrorRow>");
  qRegisterMetaType<rw::qc::ErrorRow>("rw::qc::ErrorRow");
  qRegisterMetaType<std::shared_ptr<const rw::smile::SmileSurface>>("std::shared_ptr<const rw::smile::SmileSurface>");
  wireSignals();


  // Label "project" en barre d’état
  projectLabel_ = new QLabel(this);
  projectLabel_->setObjectName("lblProjectName");
  projectLabel_->setText("No project");
  statusBar()->addPermanentWidget(projectLabel_, /*stretch*/1);

  // Raccourcis clavier pour Save / Save As
  setupProjectShortcuts();

  // initialise l’état
  setCurrentProject(QString(), false);

  setupConvergenceChart();
  setupConvergenceLegend();

  setupGreeksChart();

  // Debounce timer
  stressDebounce_ = new QTimer(this);
  stressDebounce_->setSingleShot(true);
  stressDebounce_->setInterval(300);
  connect(stressDebounce_, &QTimer::timeout, this, [this]{
    if (stressBusy_) { stressPending_ = true; return; }
    qDebug() << "[UI] stressDebounce fired → runStress()";
    onRunStress();
  });

  // Initialisations “Stress”
  initStressControls();
  setupStressChart();
  setupStressLegend();

  setupCalibrationTab_();
  wireCalibration_();


  // Valeurs par défaut raisonnables (au cas où)
  if (ui->cbType->count()==0) { ui->cbType->addItems({"Call","Put"}); }
}

MainWindow::~MainWindow() {

  // Débrancher tout ce qui part du worker (sécurité)
  if (mcWorker_) QObject::disconnect(mcWorker_, nullptr, this, nullptr);

  stopMcWorker(/*hard=*/true);                 // arrête proprement le thread si besoin

  // Détruit explicitement les vues de charts (elles possèdent leurs QChart)
  delete convChartView_;   convChartView_   = nullptr;
  delete greeksChartView_; greeksChartView_ = nullptr;
  delete stressChartView_; stressChartView_ = nullptr;

  if (calibThread_) {
    QObject::disconnect(calib_, nullptr, this, nullptr);
    calibThread_->quit();
    calibThread_->wait();
    calibThread_ = nullptr;
    calib_ = nullptr; // sera deleteLater() via finished
  }
  // Calibration charts (les QChart seront détruits par les views)
  delete chartCalibSmile_;  chartCalibSmile_ = nullptr;
  delete chartCalibResid_;  chartCalibResid_ = nullptr;
  smileChart_ = nullptr; residChart_ = nullptr;
  smileLine_ = nullptr; smileNodes_ = nullptr; residBars_ = nullptr;


  delete ui;
}

void MainWindow::wireSignals() {
  connect(ui->btnRunMC,        &QPushButton::clicked, this, &MainWindow::onRunMC);
  connect(ui->btnRunGreeks,    &QPushButton::clicked, this, &MainWindow::onRunGreeks);
  connect(ui->btnRunStress,    &QPushButton::clicked, this, &MainWindow::onRunStress);
  connect(ui->btnLoadDefaults, &QPushButton::clicked, this, &MainWindow::onLoadDefaults);
  connect(ui->btnReset,        &QPushButton::clicked, this, &MainWindow::onReset);
  // connect(ui->btnSnapBaseline, &QPushButton::clicked, this, &MainWindow::onSnapBaseline);
  connect(ui->btnSaveProject, &QPushButton::clicked, this, &MainWindow::onSaveProject);
  connect(ui->btnLoadProject, &QPushButton::clicked, this, &MainWindow::onLoadProject);
  connect(ui->btnExportPng, &QPushButton::clicked, this, &MainWindow::onExportChartsPng);

  connect(ui->chkUseSmile, &QCheckBox::toggled, this, [this](bool){
    if (mcWorker_) applySmileSettingsToWorker_();
    statusBar()->showMessage("Smile mode mis à jour.");
  });
  connect(ui->sbVegaStress, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double){
    if (mcWorker_) applySmileSettingsToWorker_();
  });

  if (ui->btnRunAllTests){connect(ui->btnRunAllTests, &QPushButton::clicked, this, &MainWindow::runAllTests);}
  


  auto tagDirtyD = [this](QDoubleSpinBox* w){
    if (!w) return;
    connect(w, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, [this]{ markProjectDirty(); });
  };
  auto tagDirtyI = [this](QSpinBox* w){
    if (!w) return;
    connect(w, qOverload<int>(&QSpinBox::valueChanged),
            this, [this]{ markProjectDirty(); });
  };

  // Doubles 
  tagDirtyD(ui->sbS0);
  tagDirtyD(ui->sbR);
  tagDirtyD(ui->sbQ);
  tagDirtyD(ui->sbSigma);
  tagDirtyD(ui->dsbTol);
  tagDirtyD(ui->sbK);
  tagDirtyD(ui->sbT);
  tagDirtyD(ui->dsbBumpRelS0);
  tagDirtyD(ui->dsbBumpAbsSigma);
  tagDirtyD(ui->dsbBumpAbsR);
  tagDirtyD(ui->dsbBumpAbsT);
  tagDirtyD(ui->sbNTarget);
  tagDirtyD(ui->sbBatch);
  
  // Ints
  tagDirtyI(ui->sbNSteps);
  tagDirtyI(ui->sbSeed);
}

void MainWindow::applySmileSettingsToWorker_() {
  if (!mcWorker_) return;

  const bool useSmile = ui->chkUseSmile && ui->chkUseSmile->isChecked();
  const double vegaPct = ui->sbVegaStress ? ui->sbVegaStress->value() : 0.0;

  // surface à pousser
  std::shared_ptr<const rw::smile::SmileSurface> surfPtr;
  if (useSmile && !smileCalibCache_.slices.empty())
    surfPtr = std::make_shared<rw::smile::SmileSurface>(smileCalibCache_);

  // tout passe en QueuedConnection
  QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, useSmile]{ w->setSmileMode(useSmile); },
      Qt::QueuedConnection);

  QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, vegaPct]{ w->setVegaStressPct(vegaPct); },
      Qt::QueuedConnection);

  QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, surfPtr]{ w->setSmileSurface(surfPtr); },
      Qt::QueuedConnection);

  qDebug() << "[UI→Worker] applySmileSettings"
           << "useSmile=" << useSmile
           << "vega%=" << vegaPct
           << "hasSurface=" << (!smileCalibCache_.slices.empty());
}


void MainWindow::setupProjectShortcuts() {
  auto* scSave = new QShortcut(QKeySequence::Save, this); // Ctrl+S
  connect(scSave, &QShortcut::activated, this, &MainWindow::onShortcutSave);

  auto* scSaveAs = new QShortcut(QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_S), this);
  connect(scSaveAs, &QShortcut::activated, this, &MainWindow::onShortcutSaveAs);
}

static void bindSliderToDoubleSpin(QSlider* sl, QDoubleSpinBox* sb,
                                   double scale, double offset = 0.0)
{
  // slider -> spinbox
  QObject::connect(sl, &QSlider::valueChanged, sb, [=](int v){
    QSignalBlocker b(sb);
    sb->setValue(offset + scale * static_cast<double>(v));
  });
  // spinbox -> slider
  QObject::connect(sb, qOverload<double>(&QDoubleSpinBox::valueChanged), sl, [=](double x){
    QSignalBlocker b(sl);
    sl->setValue(qRound((x - offset) / scale));
  });
}


/*
void MainWindow::startMcWorker() {
  stopMcWorker(); // safety

  mcThread_  = new QThread(this);
  mcWorker_  = new gui::McWorker();
  mcWorker_->moveToThread(mcThread_);

  connect(mcThread_, &QThread::finished, mcWorker_, &QObject::deleteLater);

  connect(mcWorker_, &gui::McWorker::progress, this, &MainWindow::onMcProgress);
  connect(mcWorker_, &gui::McWorker::finished, this, &MainWindow::onMcFinished);
  connect(mcWorker_, &gui::McWorker::failed,   this, &MainWindow::onMcFailed);
  connect(mcWorker_, &gui::McWorker::canceled, this, &MainWindow::onMcCanceled);

  mcThread_->start();
}
*/

void MainWindow::startMcWorker() {
  // Évite d’empiler des threads si déjà lancé
  if (mcThread_ && mcWorker_) return;

  // Pour que les signaux/queued invoke avec ces types passent entre threads
  qRegisterMetaType<std::size_t>("std::size_t");
  qRegisterMetaType<long long>("long long");
  qRegisterMetaType<rw::config::McConfig>("rw::config::McConfig");

  mcThread_  = new QThread(this);
  mcWorker_  = new gui::McWorker();             // PAS de parent → il vivra dans mcThread_
  mcWorker_->moveToThread(mcThread_);

  connect(mcThread_, &QThread::finished, mcWorker_, &QObject::deleteLater);

  // Worker -> GUI : force Qt::QueuedConnection (sécurité inter-threads)
  connect(mcWorker_, &gui::McWorker::progress, this, &MainWindow::onMcProgress, Qt::QueuedConnection);
  connect(mcWorker_, &gui::McWorker::finished, this, &MainWindow::onMcFinished, Qt::QueuedConnection);
  connect(mcWorker_, &gui::McWorker::failed,   this, &MainWindow::onMcFailed,   Qt::QueuedConnection);
  connect(mcWorker_, &gui::McWorker::canceled, this, &MainWindow::onMcCanceled, Qt::QueuedConnection);

  // (Si tu avais d’autres signaux comme greeks/stress, ajoute-les ici pareil en Queued)

  mcThread_->start();
}

/*
void MainWindow::stopMcWorkeravant() {
  if (mcWorker_) {
    QMetaObject::invokeMethod(mcWorker_, "requestStop", Qt::QueuedConnection);
  }
  if (mcThread_) {
    mcThread_->quit();
    mcThread_->wait();
    mcThread_->deleteLater();
    mcThread_ = nullptr;
    mcWorker_ = nullptr;
  }
}

*/
/*
void MainWindow::stopMcWorker(bool hard) {
  if (!mcThread_) return;

  // Demande d’arrêt côté worker
  if (mcWorker_) mcWorker_->requestStop();

  // Quitte et attend le thread
  mcThread_->quit();
  mcThread_->wait();

  // IMPORTANT : détruire synchroniquement pour que LSan ne voie pas de pending delete
  if (hard) {
    delete mcWorker_;  mcWorker_  = nullptr;
    delete mcThread_;  mcThread_  = nullptr;
  } else {
    // Chemin "normal" si l'appli continue de tourner
    mcWorker_->deleteLater(); mcWorker_ = nullptr;
    mcThread_->deleteLater(); mcThread_ = nullptr;
  }
}
*/
void MainWindow::stopMcWorker(bool hard) {
  Q_UNUSED(hard);
  if (!mcThread_) return;            // rien à faire

  if (mcWorker_) {
    // Coupe toute remontée de signaux vers 'this'
    QObject::disconnect(mcWorker_, nullptr, this, nullptr);
    // Demande d'arrêt dans le thread du worker
    QMetaObject::invokeMethod(mcWorker_, "requestStop", Qt::QueuedConnection);
  }

  mcThread_->quit();
  mcThread_->wait(3000);             // petit timeout ⇒ pas de blocage infini

  // deleteLater(worker_) est déjà connecté sur finished du thread
  mcThread_->deleteLater();
  mcThread_  = nullptr;
  mcWorker_  = nullptr;              // pointeur rendu inoffensif
}

void MainWindow::onRunMC() {
  // Lire inputs (onglet Marché & Instrument)
  const double S0    = ui->sbS0->value();
  const double r     = ui->sbR->value();
  const double q     = ui->sbQ->value();
  const double sigma = ui->sbSigma->value();

  const double K     = ui->sbK->value();
  const double T     = ui->sbT->value();
  const bool isPut   = (ui->cbType->currentText().toLower() == "put");

  rw::market::MarketData mkt(S0, r, q);
  rw::market::OptionType typ = isPut ? rw::market::OptionType::Put
                                     : rw::market::OptionType::Call;
  rw::market::Instrument inst(K, T, typ);

  // MC config (onglet Pricing)
  const std::size_t nTarget = static_cast<std::size_t>(ui->sbNTarget->value());
  const std::size_t batch   = static_cast<std::size_t>(ui->sbBatch->value());
  const double      tol     = ui->dsbTol->value();
  const std::size_t nSteps  = static_cast<std::size_t>(ui->sbNSteps->value());
  const std::uint64_t seed  = static_cast<std::uint64_t>(ui->sbSeed->value());
  
  rw::config::McConfig cfg(nTarget, batch, tol, nSteps, seed);

  // Convergence log UI
  // if (ui->cbLog && ui->cbLog->isChecked()) clearConvergenceTable();

  // clear the chart before launching
  resetConvergenceChart();
  if (ui->btnRunMC) ui->btnRunMC->setEnabled(false);

  // Lancer worker
  startMcWorker();
  applySmileSettingsToWorker_();

  // garder tout le run visible
  if (xAxis_) xAxis_->setRange(0.0, static_cast<double>(nTarget));

  // reset des bornes Y observées
  yMinSeen_ =  std::numeric_limits<double>::infinity();
  yMaxSeen_ = -std::numeric_limits<double>::infinity();


  if (bsLine_) {
    const double bs = bs_price(mkt, inst, sigma);
    bsLine_->clear();
    bsLine_->append(0.0, bs);
    bsLine_->append(static_cast<double>(nTarget), bs);
  }

  // Vide le marqueur du dernier point en début de run
  if (lastPt_) lastPt_->clear();

  
  QMetaObject::invokeMethod(
    mcWorker_,
    [w=mcWorker_, mkt, inst, cfg, sigma]() { w->runPricing(mkt, inst, cfg, sigma); },
    Qt::QueuedConnection
  );

  // Reset affichage résultats en attendant
  setPricingResults(std::nan(""), std::nan(""), std::nan(""), std::nan(""), 0, 0);
}

void MainWindow::onRunGreeks() {
  // Lire Market & Instrument
  const double S0    = ui->sbS0->value();
  const double r     = ui->sbR->value();
  const double q     = ui->sbQ->value();
  const double sigma = ui->sbSigma->value();
  const double K     = ui->sbK->value();
  const double T     = ui->sbT->value();
  const bool   isPut = (ui->cbType->currentText().toLower() == "put");

  rw::market::MarketData mkt(S0, r, q);
  rw::market::OptionType typ = isPut ? rw::market::OptionType::Put
                                     : rw::market::OptionType::Call;
  rw::market::Instrument inst(K, T, typ);

  // MC de base (N/batch/steps/seed/tol)
  const std::size_t nTarget = static_cast<std::size_t>(ui->sbNTarget->value());
  const std::size_t batch   = static_cast<std::size_t>(ui->sbBatch->value());
  const double      tol     = ui->dsbTol->value();
  const std::size_t nSteps  = static_cast<std::size_t>(ui->sbNSteps->value());
  const std::uint64_t seed  = static_cast<std::uint64_t>(ui->sbSeed->value());
  rw::config::McConfig mc(nTarget, batch, tol, nSteps, seed);

  // Méthode : BRV / PW / LRM
  QString method = "brv";
  if (ui->rbPW->isChecked())  method = "pw";
  if (ui->rbLRM->isChecked()) method = "lrm";

  // GreeksConfig (bump + anti/CRN/steps)
  rw::config::GreeksConfig gk(rw::config::Greek::Delta); // le "greek" sera choisi à l’intérieur
  gk.bump_rel_S0    = ui->dsbBumpRelS0->value();
  gk.bump_abs_sigma = ui->dsbBumpAbsSigma->value();
  gk.bump_abs_r     = ui->dsbBumpAbsR->value();
  gk.bump_abs_T     = ui->dsbBumpAbsT->value();
  gk.use_antithetic = ui->cbAnti->isChecked();
  gk.use_crn        = ui->cbCRN->isChecked();
  gk.n_steps        = nSteps;
  gk.use_pathwise   = (method == "pw"); // utile pour delta/gamma PW

  // Lancer worker dédié Greeks
  startMcWorker();
  applySmileSettingsToWorker_();

  // Connecter le signal de fin des greeks (spécifique à ce run)
  connect(mcWorker_, &gui::McWorker::greeksFinished, this,
          [this](double d,double v,double g,double r,double t,long long ms){
            setGreeksResults(d,v,g,r,t);
            updateGreeksChart(d,v,g,r,t);
            if (ui->lblElapsed) ui->lblElapsed->setText(QString::number(ms)); // optionnel
            stopMcWorker();
          });
          
  // Connecter l’échec/annulation (déjà gérés dans startMcWorker via onMcFailed/onMcCanceled)

  // Démarrer
  QMetaObject::invokeMethod(
    mcWorker_,
    [w=mcWorker_, mkt, inst, mc, gk, sigma, method]() {
      w->runGreeks(mkt, inst, mc, gk, sigma, method);
    },
    Qt::QueuedConnection
  );

  // Affichage provisoire
  setGreeksResults(0,0,0,0,0);
  statusBar()->showMessage(QString("Computing greeks (%1)...").arg(method.toUpper()));
}

void MainWindow::onLoadDefaults() {
  ui->sbS0->setValue(100.0);
  ui->sbR->setValue(0.02);
  ui->sbQ->setValue(0.0);
  ui->sbSigma->setValue(0.20);

  ui->sbK->setValue(100.0);
  ui->sbT->setValue(1.0);
  ui->cbType->setCurrentText("Call");

  ui->sbNTarget->setValue(400000);
  ui->sbBatch->setValue(100000);
  ui->dsbTol->setValue(-1.0);
  ui->sbNSteps->setValue(1);
  ui->sbSeed->setValue(42);

  // Greeks defaults
  ui->rbBRV->setChecked(true);
  ui->rbPW->setChecked(false);
  ui->rbLRM->setChecked(false);
  ui->cbCRN->setChecked(true);
  ui->cbAnti->setChecked(false);
  ui->dsbBumpRelS0->setValue(0.01);
  ui->dsbBumpAbsSigma->setValue(0.01);
  ui->dsbBumpAbsR->setValue(0.0001);
  ui->dsbBumpAbsT->setValue(1.0/365.0);
}

void MainWindow::onReset() {
  // Efface résultats pricing
  if (ui->lblPrice)  ui->lblPrice->clear();
  if (ui->lblSE)     ui->lblSE->clear();
  if (ui->lblCILow)  ui->lblCILow->clear();
  if (ui->lblCIHigh) ui->lblCIHigh->clear();
  if (ui->lblNEff)   ui->lblNEff->clear();
  if (ui->lblElapsed)ui->lblElapsed->clear();

  // Efface greeks
  setGreeksResults(0,0,0,0,0);
}

void MainWindow::setPricingResults(double price, double se, double ciLow, double ciHigh,
                                   std::size_t nEff, long long ms) {
  if (ui->lblPrice)  ui->lblPrice ->setText(std::isnan(price) ? "-" : QString::asprintf("%.6f", price));
  if (ui->lblSE)     ui->lblSE    ->setText(std::isnan(se)    ? "-" : QString::asprintf("%.9f", se));
  if (ui->lblCILow)  ui->lblCILow ->setText(std::isnan(ciLow) ? "-" : QString::asprintf("%.6f", ciLow));
  if (ui->lblCIHigh) ui->lblCIHigh->setText(std::isnan(ciHigh)? "-" : QString::asprintf("%.6f", ciHigh));
  if (ui->lblNEff)   ui->lblNEff  ->setText(nEff ? QString::number(nEff) : "-");
  if (ui->lblElapsed)ui->lblElapsed->setText(ms ? QString::number(ms) : "-");
}

void MainWindow::setGreeksResults(double d, double v, double g, double r, double t) {
  if (ui->lblDelta) ui->lblDelta->setText(QString::asprintf("%.6f", d));
  if (ui->lblVega)  ui->lblVega ->setText(QString::asprintf("%.6f", v));
  if (ui->lblGamma) ui->lblGamma->setText(QString::asprintf("%.6f", g));
  if (ui->lblRho)   ui->lblRho  ->setText(QString::asprintf("%.6f", r));
  if (ui->lblTheta) ui->lblTheta->setText(QString::asprintf("%.6f", t));
}

void MainWindow::onMcProgress(std::size_t n, double mean, double half) {
  if (meanSeries_) {
    const qreal x = static_cast<qreal>(n);
    const qreal y = static_cast<qreal>(mean);
    const qreal lo = static_cast<qreal>(mean - half);
    const qreal hi = static_cast<qreal>(mean + half);

    meanSeries_->append(x, y);
    ciUpper_->append(x, hi);
    ciLower_->append(x, lo);
    updateAxesForPoint(x, lo, hi);

    
  }
  // Marque le dernier point
  if (lastPt_) {
    lastPt_->clear();
    lastPt_->append(static_cast<qreal>(n), static_cast<qreal>(mean));
  }

  // Autoscale vertical autour de l’IC courant (petit padding)
  // Met à jour bornes Y globales (sur la CI)
  yMinSeen_ = std::min(yMinSeen_, static_cast<double>(mean - half));
  yMaxSeen_ = std::max(yMaxSeen_, static_cast<double>(mean + half));

  if (yAxis_) {
    const double span = std::max(1e-9, yMaxSeen_ - yMinSeen_);
    const double pad  = 0.10 * span; // 10% de marge visuelle
    yAxis_->setRange(yMinSeen_ - pad, yMaxSeen_ + pad);
  }
}

void MainWindow::onMcFailed(const QString& why) {
  QMessageBox::warning(this, "MC failed", why);
  if (ui->btnRunStress)    ui->btnRunStress->setEnabled(true);
  //if (ui->btnSnapBaseline) ui->btnSnapBaseline->setEnabled(true);
  unsetCursor();
  stressBusy_ = false;                                   
  if (stressPending_) { stressPending_ = false; armStressDebounce_(150); } 
  stopMcWorker();
}

void MainWindow::onMcCanceled() {
  QMessageBox::information(this, "Monte Carlo", "Run canceled.");
  if (ui->btnRunStress)    ui->btnRunStress->setEnabled(true);
  //if (ui->btnSnapBaseline) ui->btnSnapBaseline->setEnabled(true);
  unsetCursor();
  stressBusy_ = false;                                   
  if (stressPending_) { stressPending_ = false; armStressDebounce_(150); } 
  stopMcWorker();
}

void MainWindow::onMcFinished(double price, double se, double lo, double hi,
                              std::size_t nEff, long long ms) {
  setPricingResults(price, se, lo, hi, nEff, ms);
  if (ui->btnRunMC) ui->btnRunMC->setEnabled(true);
  stopMcWorker();
}

void MainWindow::onSnapBaseline() {
  const bool isPut = (ui->cbType->currentText().toLower() == "put");
  baseMkt_.emplace(
      ui->sbS0->value(),
      ui->sbR->value(),
      ui->sbQ->value()
  );

  baseInst_.emplace(
      ui->sbK->value(),
      ui->sbT->value(),
      isPut ? rw::market::OptionType::Put
            : rw::market::OptionType::Call
  );
  baseSigma_ = ui->sbSigma->value();
  baselineSet_ = true;
  statusBar()->showMessage("Baseline snapped.", 1500);
}

void MainWindow::onRunStress() {
  qDebug() << "[UI] onRunStress fire  dS%=" << (ui->sbS0Stress?ui->sbS0Stress->value():0)
         << " dSigma=" << (ui->sbSigmaStress?ui->sbSigmaStress->value():0)
         << " vega%=" << (ui->sbVegaStress?ui->sbVegaStress->value():0)
         << " useSmile=" << (ui->chkUseSmile && ui->chkUseSmile->isChecked());

  if (!baselineSet_) onSnapBaseline();   // garantit baseMkt_/baseInst_/baseSigma_

  if (stressBusy_) return;
  stressBusy_ = true;
  resetStressChart();


  // 1) Lire les chocs UI
  const double dS_rel_pct = ui->sbS0Stress    ? ui->sbS0Stress->value()    : 0.0; // en %
  const double dSigma     = ui->sbSigmaStress ? ui->sbSigmaStress->value() : 0.0; // vol pts (abs)
  const double dR_bps     = ui->sbRStress     ? ui->sbRStress->value()     : 0.0; // bps
  const double dT_days    = ui->sbTStress     ? ui->sbTStress->value()     : 0.0; // jours

  // 2) Conversions d’unités
  const double dS_rel   = dS_rel_pct / 100.0;
  const double dR_abs   = dR_bps     / 10000.0;
  const double dT_years = dT_days    / 365.0;

  // (optionnel) garde-fous simples côté UI
  if (baseMkt_) {
    const double S0s = baseMkt_->S0 * (1.0 + dS_rel);
    if (S0s <= 0.0) {
      QMessageBox::warning(this, "Stress", "ΔS trop négatif (S0 stressé ≤ 0).");
      return;
    }
  }
  if (baseSigma_ + dSigma < 0.0) {
    QMessageBox::warning(this, "Stress", "σ stressée < 0.");
    return;
  }

  // 3) McConfig (réutilise l’onglet Pricing)
  const std::size_t nTarget = static_cast<std::size_t>(ui->sbNTarget->value());
  const std::size_t batch   = static_cast<std::size_t>(ui->sbBatch->value());
  const double      tol     = ui->dsbTol->value();
  const std::size_t nSteps  = static_cast<std::size_t>(ui->sbNSteps->value());
  const std::uint64_t seed  = static_cast<std::uint64_t>(ui->sbSeed->value());
  rw::config::McConfig cfg(nTarget, batch, tol, nSteps, seed);

  // 4) UI: état “busy”
  if (ui->btnRunStress)   ui->btnRunStress->setEnabled(false);
  // if (ui->btnSnapBaseline)  ui->btnSnapBaseline->setEnabled(false);
  setCursor(Qt::BusyCursor);

  
  startMcWorker();
  applySmileSettingsToWorker_();


  // 5) Approx Taylor si grecs déjà présents dans l’UI
  auto readLbl = [](QLabel* L, double& out) {
    if (!L) return false;
    bool ok=false; const double v = L->text().toDouble(&ok);
    if (ok) out=v; return ok;
  };

  // on les calcule avant, et on les réutilisera dans le handler
  double d_contrib=0.0, g_contrib=0.0, v_contrib=0.0, r_contrib=0.0;

  if (baseMkt_) {
    double Delta=0, Gamma=0, Vega=0, Rho=0;
    const bool hasD = readLbl(ui->lblDelta, Delta);
    const bool hasG = readLbl(ui->lblGamma, Gamma);
    const bool hasV = readLbl(ui->lblVega,  Vega);
    const bool hasR = readLbl(ui->lblRho,   Rho);

    if (hasD || hasG || hasV || hasR) {
      const double dS_abs = baseMkt_->S0 * dS_rel;
      if (hasD) d_contrib = Delta * dS_abs;
      if (hasG) g_contrib = 0.5  * Gamma * dS_abs * dS_abs;
      if (hasV) v_contrib = Vega  * dSigma;
      if (hasR) r_contrib = Rho   * dR_abs;

      const double pnlTaylor = d_contrib + g_contrib + v_contrib + r_contrib;
      if (ui->lblPLTaylor) ui->lblPLTaylor->setText(QString::asprintf("%.6f", pnlTaylor));

      // on met à jour tout de suite le graphe côté Taylor (barres à gauche)
      updateStressChart(/*pnlMC=*/0.0, d_contrib, g_contrib, v_contrib, r_contrib);
    } else {
      if (ui->lblPLTaylor) ui->lblPLTaylor->setText("-");
      // on vide la partie Taylor du graphe (affiche 0)
      updateStressChart(/*pnlMC=*/0.0, 0.0, 0.0, 0.0, 0.0);
      statusBar()->showMessage("Calcule d’abord les Greeks pour voir l’approx Taylor.", 2500);
    }
  }

  // 6) Handler de fin (unique)
  connect(mcWorker_, &gui::McWorker::stressFinished, this,
          //  on CAPTURE les contributions par valeur pour pouvoir dessiner la barre MC en face
          [this, d_contrib, g_contrib, v_contrib, r_contrib]
          (double baseP, double stressP, double pnlMC, long long ms) {

            if (ui->lblPL)        ui->lblPL->setText(QString::asprintf("%.6f", pnlMC));
            if (ui->lblPrice)     ui->lblPrice->setText(QString::asprintf("%.6f", stressP));
            if (ui->lblElapsed)   ui->lblElapsed->setText(QString::number(ms));

            // Erreur vs Taylor si présent
            if (ui->lblPLTaylor && !ui->lblPLTaylor->text().isEmpty()) {
              bool ok=false; const double t = ui->lblPLTaylor->text().toDouble(&ok);
              if (ok && ui->lblPLError) ui->lblPLError->setText(QString::asprintf("%.6f", pnlMC - t));
            }

            //  met à jour le graphe côté MC (barre à droite)
            updateStressChart(/*pnlMC=*/pnlMC, d_contrib, g_contrib, v_contrib, r_contrib);

            // UI: fin “busy”
            if (ui->btnRunStress)    ui->btnRunStress->setEnabled(true);
            // if (ui->btnSnapBaseline) ui->btnSnapBaseline->setEnabled(true);
            unsetCursor();

            statusBar()->showMessage("Stress done.", 1500);
            stressBusy_ = false; 
            if (stressPending_) { stressPending_ = false; armStressDebounce_(150); } 
            stopMcWorker();
          },
          Qt::UniqueConnection);

  // 7) Lancer le run stress (déréférencer les optionnels)
  QMetaObject::invokeMethod(
    mcWorker_,
    [w=mcWorker_, this, cfg, dS_rel, dSigma, dR_abs, dT_years]() {
      w->runStress(*baseMkt_, *baseInst_, cfg, baseSigma_, dS_rel, dSigma, dR_abs, dT_years);
    },
    Qt::QueuedConnection
  );

  statusBar()->showMessage("Running stress…");
}

// ======================= Convergence =======================
void MainWindow::setupConvergenceChart() {
  if (convChartView_) return;

  using namespace QtCharts;

  // Le chart sera le parent de toutes les séries/axes
  auto* chart = new QChart();
  chart->legend()->setVisible(false);
  chart->setTitle("MC Convergence");

  // Axes (parent = chart)
  xAxis_ = new QValueAxis(chart);
  yAxis_ = new QValueAxis(chart);
  xAxis_->setTitleText("N");
  yAxis_->setTitleText("Price");
  xAxis_->setLabelFormat("%.0f");
  yAxis_->setLabelFormat("%.6f");
  chart->addAxis(xAxis_, Qt::AlignBottom);
  chart->addAxis(yAxis_, Qt::AlignLeft);

  // Séries (parent = chart)
  meanSeries_ = new QLineSeries(chart);   meanSeries_->setName("MC mean");
  ciUpper_    = new QLineSeries(chart);
  ciLower_    = new QLineSeries(chart);
  ciBand_     = new QAreaSeries(ciUpper_, ciLower_); 
  ciBand_->setParent(chart);              // QAreaSeries n’a pas de parent dans le ctor

  bsLine_     = new QLineSeries(chart);   bsLine_->setName("BS ref");
  lastPt_     = new QScatterSeries(chart); 
  lastPt_->setName("Last batch");

  // Styles
  meanSeries_->setPen(QPen(C_MEAN, 2));

  ciUpper_->setPen(Qt::NoPen);
  ciLower_->setPen(Qt::NoPen);
  ciBand_->setPen(Qt::NoPen);
  ciBand_->setBrush(QColor(C_CI.red(), C_CI.green(), C_CI.blue(), A_CI));

  bsLine_->setPen(QPen(C_BS, 1.8, Qt::SolidLine));

  lastPt_->setMarkerShape(QScatterSeries::MarkerShapeCircle);
  lastPt_->setMarkerSize(9.0);
  lastPt_->setColor(C_LAST);
  lastPt_->setBorderColor(C_LAST);

  // Ajout + attache axes
  chart->addSeries(ciBand_);
  chart->addSeries(meanSeries_);
  chart->addSeries(bsLine_);
  chart->addSeries(lastPt_);

  ciBand_->attachAxis(xAxis_);     ciBand_->attachAxis(yAxis_);
  meanSeries_->attachAxis(xAxis_); meanSeries_->attachAxis(yAxis_);
  bsLine_->attachAxis(xAxis_);     bsLine_->attachAxis(yAxis_);
  lastPt_->attachAxis(xAxis_);     lastPt_->attachAxis(yAxis_);

  // ChartView parentée au container
  convChartView_ = new QChartView(chart, ui->chartContainer);
  convChartView_->setRenderHint(QPainter::Antialiasing);

  auto* lay = new QVBoxLayout(ui->chartContainer);
  lay->setContentsMargins(0,0,0,0);
  lay->addWidget(convChartView_);
}

void MainWindow::resetConvergenceChart() {
  setupConvergenceChart();
  meanSeries_->clear();
  ciUpper_->clear();
  ciLower_->clear();
  if (ciBand_) { ciBand_->setUpperSeries(ciUpper_); ciBand_->setLowerSeries(ciLower_); }
  if (lastPt_) lastPt_->clear();
}

void MainWindow::updateAxesForPoint(qreal x, qreal yLow, qreal yHigh) {
  if (!xAxis_ || !yAxis_) return;

  if (meanSeries_->count() == 1) {
    // first point: set tight range
    xAxis_->setRange(0, std::max<qreal>(1.0, x));
    yAxis_->setRange(std::min(yLow, yHigh), std::max(yLow, yHigh));
  } else {
    if (x > xAxis_->max())  xAxis_->setMax(x);
    if (yLow  < yAxis_->min()) yAxis_->setMin(yLow);
    if (yHigh > yAxis_->max()) yAxis_->setMax(yHigh);
  }
}

static QWidget* makeLineSwatch(const QColor& c, int w=26, int h=12,
                               Qt::PenStyle style=Qt::SolidLine) {
  auto* L = new QLabel; L->setFixedSize(w,h);
  QPixmap pm(w,h); pm.fill(Qt::transparent);
  QPainter p(&pm); p.setRenderHint(QPainter::Antialiasing, true);
  QPen pen(c, 2, style); p.setPen(pen);
  p.drawLine(2, h/2, w-2, h/2);
  L->setPixmap(pm); return L;
}

static QWidget* makeAreaSwatch(const QColor& c, int w=26, int h=12) {
  auto* L = new QLabel; L->setFixedSize(w,h);
  QPixmap pm(w,h); pm.fill(Qt::transparent);
  QPainter p(&pm); p.setRenderHint(QPainter::Antialiasing, true);
  p.fillRect(QRect(2,2,w-4,h-4), c);
  L->setPixmap(pm); return L;
}

static QWidget* makeDotSwatch(const QColor& c, int size=10) {
  auto* L = new QLabel; L->setFixedSize(size+4,size+4);
  QPixmap pm(size+4,size+4); pm.fill(Qt::transparent);
  QPainter p(&pm); p.setRenderHint(QPainter::Antialiasing, true);
  p.setPen(Qt::NoPen); p.setBrush(c);
  p.drawEllipse(QPoint((size+4)/2,(size+4)/2), size/2, size/2);
  L->setPixmap(pm); return L;
}

void MainWindow::setupConvergenceLegend() {
  if (!ui->legendLayoutPricing) return;  // le QVBoxLayout dans ton QGroupBox

  // Nettoyage du contenu actuel du layout (récursif)
  auto clearLayout = [](QLayout* layout, auto&& clearRef) -> void {
    while (QLayoutItem* it = layout->takeAt(0)) {
      if (auto* w = it->widget()) { w->deleteLater(); }
      if (auto* l = it->layout()) { clearRef(l, clearRef); delete l; }
      delete it;
    }
  };
  clearLayout(ui->legendLayoutPricing, clearLayout);

  // Ajoute une ligne (swatch + label) au VLayout via un HLayout
  auto addRow = [&](QWidget* sw, const QString& txt){
    auto* row = new QHBoxLayout();
    row->setContentsMargins(0,0,0,0);
    row->setSpacing(8);
    auto* label = new QLabel(txt);
    row->addWidget(sw,   0, Qt::AlignLeft | Qt::AlignVCenter);
    row->addWidget(label,0, Qt::AlignLeft | Qt::AlignVCenter);
    ui->legendLayoutPricing->addLayout(row);
  };

  // Couleurs/styles cohérents avec le chart
  addRow(makeLineSwatch(C_MEAN),                                     tr("MC mean"));
  addRow(makeAreaSwatch(QColor(C_CI.red(), C_CI.green(), C_CI.blue(), A_CI)), tr("95% CI band"));
  addRow(makeLineSwatch(C_BS, 26, 12, Qt::SolidLine),                tr("BS ref"));
  addRow(makeDotSwatch (C_LAST, 10),                                 tr("Last batch"));

  ui->legendLayoutPricing->addStretch(); // pousse tout en haut du groupbox
}

QJsonObject MainWindow::makeProjectJson() const {
  // market
  QJsonObject market{
    {"S0",    ui->sbS0->value()},
    {"r",     ui->sbR->value()},
    {"q",     ui->sbQ->value()},
    {"sigma", ui->sbSigma->value()}
  };

  // instrument
  const bool isPut = (ui->cbType->currentText().toLower() == "put");
  QJsonObject instrument{
    {"type", isPut ? "put" : "call"},
    {"K", ui->sbK->value()},
    {"T", ui->sbT->value()}
  };

  // MC config (+ flags anti/CRN de l’onglet Greeks)
  QJsonObject mc{
    {"seed",            static_cast<double>(ui->sbSeed->value())},
    {"n_paths_target",  static_cast<double>(ui->sbNTarget->value())},
    {"batch_size",      static_cast<double>(ui->sbBatch->value())},
    {"tolerance",       ui->dsbTol->value()},
    {"n_steps",         static_cast<double>(ui->sbNSteps->value())},
    {"use_antithetic",  ui->cbAnti ? ui->cbAnti->isChecked() : false},
    {"use_crn",         ui->cbCRN  ? ui->cbCRN ->isChecked() : true}
  };

  // UI info (onglet courant)
  QString lastTab = "Pricing";
  if (ui->tabWidget)
    lastTab = ui->tabWidget->tabText(ui->tabWidget->currentIndex());

  QJsonObject uiobj{ {"last_tab", lastTab} };

  QJsonObject root;
  root["market"]     = market;
  root["instrument"] = instrument;
  root["mc"]         = mc;
  root["ui"]         = uiobj;

    // ============= Greeks =============
  // config greeks (méthode + bumps)
  QString gMethod = "brv";
  if (ui->rbPW && ui->rbPW->isChecked()) gMethod = "pw";
  if (ui->rbLRM && ui->rbLRM->isChecked()) gMethod = "lrm";

  QJsonObject greeks_cfg{
    {"method",         gMethod},
    {"bump_rel_S0",    ui->dsbBumpRelS0 ? ui->dsbBumpRelS0->value() : 0.0},
    {"bump_abs_sigma", ui->dsbBumpAbsSigma ? ui->dsbBumpAbsSigma->value() : 0.0},
    {"bump_abs_r",     ui->dsbBumpAbsR ? ui->dsbBumpAbsR->value() : 0.0},
    {"bump_abs_T",     ui->dsbBumpAbsT ? ui->dsbBumpAbsT->value() : 0.0},
    // flags déjà dans "mc", on les duplique ici pour tracer l’historique exact
    {"use_antithetic", ui->cbAnti ? ui->cbAnti->isChecked() : false},
    {"use_crn",        ui->cbCRN  ? ui->cbCRN ->isChecked() : true},
    {"n_steps",        static_cast<double>(ui->sbNSteps ? ui->sbNSteps->value() : 1)}
  };

  // résultats greeks (labels)
  QJsonObject greeks_res{
    {"delta", labelToJson(ui->lblDelta)},
    {"vega",  labelToJson(ui->lblVega)},
    {"gamma", labelToJson(ui->lblGamma)},
    {"rho",   labelToJson(ui->lblRho)},
    {"theta", labelToJson(ui->lblTheta)}
  };

  QJsonObject greeks{
    {"config",  greeks_cfg},
    {"result",  greeks_res}
  };
  root["greeks"] = greeks;

  // ============= Stress =============
  // baseline (si disponible)
  QJsonObject baseline;
  if (baselineSet_ && baseMkt_ && baseInst_) {
    baseline = QJsonObject{
      {"market",      QJsonObject{{"S0", baseMkt_->S0}, {"r", baseMkt_->r}, {"q", baseMkt_->q}}},
      {"instrument",  QJsonObject{{"type", (baseInst_->type==rw::market::OptionType::Put ? "put" : "call")},
                                  {"K", baseInst_->K}, {"T", baseInst_->T}}},
      {"sigma",       baseSigma_}
    };
  }

  // chocs (UI) + conversions
  const double dS_pct  = ui->sbS0Stress    ? ui->sbS0Stress->value()    : 0.0;
  const double dSigma  = ui->sbSigmaStress ? ui->sbSigmaStress->value() : 0.0;
  const double dR_bps  = ui->sbRStress     ? ui->sbRStress->value()     : 0.0;
  const double dT_days = ui->sbTStress     ? ui->sbTStress->value()     : 0.0;

  QJsonObject stress_cfg{
    {"dS_pct",  dS_pct},
    {"dSigma",  dSigma},
    {"dR_bps",  dR_bps},
    {"dT_days", dT_days},
    // conversions qui ont servi au run
    {"derived", QJsonObject{
       {"dS_rel",   dS_pct/100.0},
       {"dR_abs",   dR_bps/10000.0},
       {"dT_years", dT_days/365.0}
    }},
    // rappel MC pour traçabilité
    {"seed",        static_cast<double>(ui->sbSeed ? ui->sbSeed->value() : 0)},
    {"n_steps",     static_cast<double>(ui->sbNSteps ? ui->sbNSteps->value() : 1)},
    {"use_crn",     ui->cbCRN ? ui->cbCRN->isChecked() : true}
  };

  // résultats stress (labels affichés)
  QJsonObject stress_res{
    {"pnl_mc",      labelToJson(ui->lblPL)},
    {"pnl_taylor",  labelToJson(ui->lblPLTaylor)},
    {"pnl_error",   labelToJson(ui->lblPLError)},
    {"stressed_px", labelToJson(ui->lblPriceStressed)},
    {"elapsed_ms",  labelToJson(ui->lblElapsed)}
  };

  QJsonObject stress{
    {"baseline", baseline},     // peut être vide/null
    {"config",   stress_cfg},
    {"result",   stress_res}
  };
  root["stress"] = stress;

  // ============= Smile / Surface =============
  QJsonObject smile;
  const bool useSmile = (ui->chkUseSmile ? ui->chkUseSmile->isChecked() : false);
  const double vegaPct = (ui->sbVegaStress ? ui->sbVegaStress->value() : 0.0);

  // si tu as l’underlying côté CalibWorker, ajoute-le aussi
  smile["enabled"] = useSmile;
  smile["vega_stress_pct"] = vegaPct;

  // surface sérialisée si dispo
  if (!smileCalibCache_.slices.empty()) {
    smile["slices"] = serializeSmileSlices_();
    // Optionnel: bornes, compte de points, time-stamp
    smile["n_slices"] = static_cast<int>(smileCalibCache_.slices.size());
  } else {
    smile["slices"] = QJsonArray(); // vide
  }

  root["smile"] = smile;
  return root;
}

void MainWindow::loadProjectJson(const QJsonObject& root) {
  // market
  if (auto m = root["market"].toObject(); !m.isEmpty()) {
    if (ui->sbS0)    ui->sbS0->setValue(m["S0"].toDouble(ui->sbS0->value()));
    if (ui->sbR)     ui->sbR ->setValue(m["r"].toDouble(ui->sbR->value()));
    if (ui->sbQ)     ui->sbQ ->setValue(m["q"].toDouble(ui->sbQ->value()));
    if (ui->sbSigma) ui->sbSigma->setValue(m["sigma"].toDouble(ui->sbSigma->value()));
  }

  // instrument
  if (auto ins = root["instrument"].toObject(); !ins.isEmpty()) {
    if (ui->sbK) ui->sbK->setValue(ins["K"].toDouble(ui->sbK->value()));
    if (ui->sbT) ui->sbT->setValue(ins["T"].toDouble(ui->sbT->value()));
    if (ui->cbType) {
      const QString t = ins["type"].toString().toLower();
      ui->cbType->setCurrentText(t == "put" ? "Put" : "Call");
    }
  }

  // MC
  if (auto mc = root["mc"].toObject(); !mc.isEmpty()) {
    if (ui->sbSeed)    ui->sbSeed   ->setValue(static_cast<int>(mc["seed"].toDouble(ui->sbSeed->value())));
    if (ui->sbNTarget) ui->sbNTarget->setValue(static_cast<int>(mc["n_paths_target"].toDouble(ui->sbNTarget->value())));
    if (ui->sbBatch)   ui->sbBatch  ->setValue(static_cast<int>(mc["batch_size"].toDouble(ui->sbBatch->value())));
    if (ui->dsbTol)    ui->dsbTol   ->setValue(mc["tolerance"].toDouble(ui->dsbTol->value()));
    if (ui->sbNSteps)  ui->sbNSteps ->setValue(static_cast<int>(mc["n_steps"].toDouble(ui->sbNSteps->value())));
    if (ui->cbAnti)    ui->cbAnti   ->setChecked(mc["use_antithetic"].toBool(ui->cbAnti->isChecked()));
    if (ui->cbCRN)     ui->cbCRN    ->setChecked(mc["use_crn"].toBool(ui->cbCRN->isChecked()));
  }
    // ============= Greeks =============
  if (auto g = root["greeks"].toObject(); !g.isEmpty()) {
    if (auto cfg = g["config"].toObject(); !cfg.isEmpty()) {
      const QString m = cfg["method"].toString();
      if (ui->rbBRV && ui->rbPW && ui->rbLRM) {
        ui->rbBRV->setChecked(m=="brv" || m.isEmpty());
        ui->rbPW ->setChecked(m=="pw");
        ui->rbLRM->setChecked(m=="lrm");
      }
      if (ui->dsbBumpRelS0)    ui->dsbBumpRelS0->setValue(cfg["bump_rel_S0"].toDouble(ui->dsbBumpRelS0->value()));
      if (ui->dsbBumpAbsSigma) ui->dsbBumpAbsSigma->setValue(cfg["bump_abs_sigma"].toDouble(ui->dsbBumpAbsSigma->value()));
      if (ui->dsbBumpAbsR)     ui->dsbBumpAbsR->setValue(cfg["bump_abs_r"].toDouble(ui->dsbBumpAbsR->value()));
      if (ui->dsbBumpAbsT)     ui->dsbBumpAbsT->setValue(cfg["bump_abs_T"].toDouble(ui->dsbBumpAbsT->value()));
      if (ui->cbAnti)          ui->cbAnti->setChecked(cfg["use_antithetic"].toBool(ui->cbAnti->isChecked()));
      if (ui->cbCRN)           ui->cbCRN ->setChecked(cfg["use_crn"].toBool(ui->cbCRN->isChecked()));
      if (ui->sbNSteps)        ui->sbNSteps->setValue(static_cast<int>(cfg["n_steps"].toDouble(ui->sbNSteps->value())));
    }
    if (auto res = g["result"].toObject(); !res.isEmpty()) {
      setLabelFromJson(ui->lblDelta, res, "delta");
      setLabelFromJson(ui->lblVega,  res, "vega");
      setLabelFromJson(ui->lblGamma, res, "gamma");
      setLabelFromJson(ui->lblRho,   res, "rho");
      setLabelFromJson(ui->lblTheta, res, "theta");
      // (optionnel) mettez à jour votre graphe de greeks si vous en avez un
      // updateGreeksChart(...);
    }
  }

  // ============= Stress =============
  if (auto s = root["stress"].toObject(); !s.isEmpty()) {
    // baseline (si présent)
    if (auto b = s["baseline"].toObject(); !b.isEmpty()) {
      auto bm = b["market"].toObject();
      auto bi = b["instrument"].toObject();
      if (!bm.isEmpty() && !bi.isEmpty() && b.contains("sigma")) {
        const double S0b = bm["S0"].toDouble(100.0);
        const double rb  = bm["r"].toDouble(0.02);
        const double qb  = bm["q"].toDouble(0.0);
        const double Kb  = bi["K"].toDouble(100.0);
        const double Tb  = bi["T"].toDouble(1.0);
        const bool isPutB = (bi["type"].toString().toLower()=="put");
        baseMkt_.emplace(S0b, rb, qb);
        baseInst_.emplace(Kb, Tb, isPutB ? rw::market::OptionType::Put
                                         : rw::market::OptionType::Call);
        baseSigma_   = b["sigma"].toDouble(0.20);
        baselineSet_ = true;
      }
    }

    // config (chocs)
    if (auto cfg = s["config"].toObject(); !cfg.isEmpty()) {
      if (ui->sbS0Stress)    ui->sbS0Stress   ->setValue(cfg["dS_pct"].toDouble(ui->sbS0Stress->value()));
      if (ui->sbSigmaStress) ui->sbSigmaStress->setValue(cfg["dSigma"].toDouble(ui->sbSigmaStress->value()));
      if (ui->sbRStress)     ui->sbRStress    ->setValue(cfg["dR_bps"].toDouble(ui->sbRStress->value()));
      if (ui->sbTStress)     ui->sbTStress    ->setValue(cfg["dT_days"].toDouble(ui->sbTStress->value()));
      // si vous avez des sliders liés, ils se resynchroniseront via vos bindings
    }

    // résultats (labels)
    if (auto res = s["result"].toObject(); !res.isEmpty()) {
      setLabelFromJson(ui->lblPL,       res, "pnl_mc");
      setLabelFromJson(ui->lblPLTaylor, res, "pnl_taylor");
      setLabelFromJson(ui->lblPLError,  res, "pnl_error");
      setLabelFromJson(ui->lblPriceStressed,    res, "stressed_px");
      setLabelFromJson(ui->lblElapsed,  res, "elapsed_ms", /*prec*/0);
    }
  }

  // UI
  if (auto u = root["ui"].toObject(); !u.isEmpty()) {
    const QString wanted = u["last_tab"].toString();
    if (ui->tabWidget && !wanted.isEmpty()) {
      for (int i=0;i<ui->tabWidget->count();++i) {
        if (ui->tabWidget->tabText(i).compare(wanted, Qt::CaseInsensitive)==0) {
          ui->tabWidget->setCurrentIndex(i);
          break;
        }
      }
    }
  }
}

void MainWindow::onSaveProject() {
  const QString dir = projectsDir();
  QDir().mkpath(dir);  // au cas où

  const QString suggested = dir + "/project_" +
      QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss") + ".json";

  const QString fn = QFileDialog::getSaveFileName(
      this, tr("Save project"), suggested, tr("JSON (*.json)"));
  if (fn.isEmpty()) return;

  QFile f(fn);
  if (!f.open(QIODevice::WriteOnly)) {
    QMessageBox::warning(this, tr("Save project"), tr("Cannot open file for writing."));
    return;
  }
  QJsonDocument doc(makeProjectJson());
  f.write(doc.toJson(QJsonDocument::Indented));
  f.close();
  setCurrentProject(fn, /*dirty*/false);
  statusBar()->showMessage(tr("Project saved to %1").arg(QDir::toNativeSeparators(fn)), 2000);
}

void MainWindow::onLoadProject() {
  const QString dir = projectsDir();
  const QString fn = QFileDialog::getOpenFileName(
      this, tr("Load project"), dir, tr("JSON (*.json)"));
  if (fn.isEmpty()) return;

  QFile f(fn);
  if (!f.open(QIODevice::ReadOnly)) {
    QMessageBox::warning(this, tr("Load project"), tr("Cannot open file for reading."));
    return;
  }
  const QByteArray bytes = f.readAll(); f.close();
  const QJsonDocument doc = QJsonDocument::fromJson(bytes);
  if (!doc.isObject()) {
    QMessageBox::warning(this, tr("Load project"), tr("Invalid JSON file."));
    return;
  }
  loadProjectJson(doc.object());
  setCurrentProject(fn, /*dirty*/false);
  qDebug() << "[UI] project loaded -> slices =" << smileCalibCache_.slices.size();
  if (ui->chkUseSmile->isChecked())
    applySmileSettingsToWorker_();
  statusBar()->showMessage(tr("Project loaded from %1").arg(QDir::toNativeSeparators(fn)), 2000);
}

QString MainWindow::projectsDir() const {
  // Point de départ = dossier de l'exécutable (…/build/bin)
  QDir exeDir(QCoreApplication::applicationDirPath());

  // Candidats (selon où tu lances le binaire)
  const QStringList candidates = {
    exeDir.absoluteFilePath("../data/projects"),
    exeDir.absoluteFilePath("../../data/projects"),
    QDir::current().absoluteFilePath("data/projects")
  };

  for (const QString& c : candidates) {
    QDir dir(c);
    if (dir.exists())
      return QDir::cleanPath(c);
  }

  // Si aucun n’existe, on crée le premier
  const QString target = QDir::cleanPath(candidates.front());
  QDir().mkpath(target);
  return target;
}

QString MainWindow::projectDisplayName() const {
  if (currentProjectPath_.isEmpty()) return "No project";
  return QFileInfo(currentProjectPath_).fileName()
         + (projectDirty_ ? " *" : "");
}

void MainWindow::setCurrentProject(const QString& path, bool dirty) {
  currentProjectPath_ = path;
  projectDirty_ = dirty;
  if (projectLabel_) projectLabel_->setText(projectDisplayName());

  // Met à jour le titre de la fenêtre
  QString base = "RiskWorkbench";
  if (!currentProjectPath_.isEmpty())
    base += " — " + projectDisplayName();
  setWindowTitle(base);
}

void MainWindow::markProjectDirty() {
  if (!projectDirty_) {
    projectDirty_ = true;
    if (projectLabel_) projectLabel_->setText(projectDisplayName());
    // mets aussi à jour le titre
    QString base = windowTitle();
    setWindowTitle("RiskWorkbench — " + projectDisplayName());
  }
}

bool MainWindow::writeProjectTo(const QString& path, QString* errMsg) const {
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    if (errMsg) *errMsg = f.errorString();
    return false;
  }
  QJsonDocument doc(makeProjectJson());
  f.write(doc.toJson(QJsonDocument::Indented));
  f.close();
  return true;
}


void MainWindow::onShortcutSave() {
  // Si on a déjà un projet ouvert => on écrase sans demander
  if (!currentProjectPath_.isEmpty()) {
    QString err;
    if (!writeProjectTo(currentProjectPath_, &err)) {
      QMessageBox::warning(this, tr("Save project"), err);
      return;
    }
    setCurrentProject(currentProjectPath_, /*dirty*/false);
    statusBar()->showMessage(tr("Project saved."), 1500);
    return;
  }

  // Pas de projet courant => “Save As…”
  onShortcutSaveAs();
}

void MainWindow::onShortcutSaveAs() {
  // Réutilise ton “Save As…” existant
  onSaveProject();
}

using namespace QtCharts;

// ========================= Greeks =========================
void MainWindow::setupGreeksChart() {
  if (greeksChartView_) return;
  if (!ui->greeksChartContainer) return;

  using namespace QtCharts;

  greeksChart_ = new QChart();
  greeksChart_->setTitle("Greeks");
  greeksChart_->legend()->setVisible(false);

  // Axes (parent = chart)
  QStringList cats = {"Delta (€/€)", "Vega (€/σ)", "Gamma (€/€²)", "Rho (€/r)", "Theta (€/an)"};
  gAxisX_ = new QBarCategoryAxis(greeksChart_); gAxisX_->append(cats);
  gAxisY_ = new QValueAxis(greeksChart_);
  gAxisY_->setTitleText("Value");
  gAxisY_->setLabelFormat("%.6f");

  // Série barres (parentages)
  greeksSeries_ = new QBarSeries(greeksChart_);
  greeksSet_    = new QBarSet("Greeks", greeksSeries_);
  greeksSet_->append({0.0,0.0,0.0,0.0,0.0});
  greeksSeries_->append(greeksSet_);

  // Style
  greeksSet_->setBorderColor(Qt::transparent);
  greeksSeries_->setBarWidth(0.6);

  greeksChart_->addSeries(greeksSeries_);
  greeksChart_->addAxis(gAxisX_, Qt::AlignBottom);
  greeksChart_->addAxis(gAxisY_, Qt::AlignLeft);
  greeksSeries_->attachAxis(gAxisX_);
  greeksSeries_->attachAxis(gAxisY_);

  // Ligne zéro (parent = chart)
  auto* zeroLine = new QLineSeries(greeksChart_);
  zeroLine->setPen(QPen(Qt::gray, 1, Qt::DashLine));
  zeroLine->append(0.0, 0.0);
  zeroLine->append(4.0, 0.0);
  greeksChart_->addSeries(zeroLine);
  zeroLine->attachAxis(gAxisX_);
  zeroLine->attachAxis(gAxisY_);

  // Vue parentée au container
  greeksChartView_ = new QChartView(greeksChart_, ui->greeksChartContainer);
  greeksChartView_->setRenderHint(QPainter::Antialiasing);

  auto* lay = new QVBoxLayout(ui->greeksChartContainer);
  lay->setContentsMargins(0,0,0,0);
  lay->addWidget(greeksChartView_);

  resetGreeksChart();
}

// Remet le chart des greeks à zéro
void MainWindow::resetGreeksChart() {
  if (!greeksSet_) return;
  greeksMinSeen_ = 0.0;
  greeksMaxSeen_ = 0.0;

  // Vide le set puis pousse 5 zéros
  const int n = greeksSet_->count();
  if (n > 0) greeksSet_->remove(0, n);
  *greeksSet_ << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;

  // Option : recaler l’axe Y
  if (gAxisY_) gAxisY_->setRange(-1.0, 1.0);
}


// Met à jour les barres (Delta, Vega, Gamma, Rho, Theta)
void MainWindow::updateGreeksChart(double d, double v, double g, double r, double t) {
  if (!greeksSet_) return;

  // Méthode robuste : on remplace tout le contenu
  const int n = greeksSet_->count();
  if (n > 0) greeksSet_->remove(0, n);
  *greeksSet_ << d << v << g << r << t;

  // Autoscale vertical simple avec un peu de marge
  greeksMinSeen_ = std::min({greeksMinSeen_, d, v, g, r, t, 0.0});
  greeksMaxSeen_ = std::max({greeksMaxSeen_, d, v, g, r, t, 0.0});
  const double pad = 0.05 * std::max(std::abs(greeksMinSeen_), std::abs(greeksMaxSeen_));
  if (gAxisY_) gAxisY_->setRange(greeksMinSeen_ - pad, greeksMaxSeen_ + pad);

}

void MainWindow::initStressControls() {
  // Ranges déjà posés dans le Designer – on relie et on mappe.

  // ΔS (%) : direct
  connect(ui->slS0Stress, &QSlider::valueChanged, this, [this](int v){
    if (ui->sbS0Stress && ui->sbS0Stress->value()!=v) ui->sbS0Stress->setValue(v);
    armStressDebounce_();
  });
  connect(ui->sbS0Stress, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v){
    if (ui->slS0Stress && ui->slS0Stress->value()!=static_cast<int>(std::round(v)))
      ui->slS0Stress->setValue(static_cast<int>(std::round(v)));
    armStressDebounce_();
  });

  // Δσ (abs) : slider −100..100 -> step 0.001
  auto sigmaFromSlider = [](int s){ return s * 0.001; };
  auto sliderFromSigma = [](double x){ return static_cast<int>(std::round(x / 0.001)); };
  connect(ui->slSigmaStress, &QSlider::valueChanged, this, [this, sigmaFromSlider](int s){
    const double v = sigmaFromSlider(s);
    if (ui->sbSigmaStress && std::abs(ui->sbSigmaStress->value()-v) > 1e-9) ui->sbSigmaStress->setValue(v);
    armStressDebounce_();
  });
  connect(ui->sbSigmaStress, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this, sliderFromSigma](double v){
    const int s = sliderFromSigma(v);
    if (ui->slSigmaStress && ui->slSigmaStress->value()!=s) ui->slSigmaStress->setValue(s);
    armStressDebounce_();
  });

  // Δr (bps) : direct
  connect(ui->slRStress, &QSlider::valueChanged, this, [this](int v){
    if (ui->sbRStress && ui->sbRStress->value()!=v) ui->sbRStress->setValue(v);
    armStressDebounce_();
  });
  connect(ui->sbRStress, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v){
    if (ui->slRStress && ui->slRStress->value()!=static_cast<int>(v))
      ui->slRStress->setValue(static_cast<int>(v));
    armStressDebounce_();
  });

  // ΔT (jours) : direct
  connect(ui->slTStress, &QSlider::valueChanged, this, [this](int v){
    if (ui->sbTStress && ui->sbTStress->value()!=v) ui->sbTStress->setValue(v);
    armStressDebounce_();
  });
  connect(ui->sbTStress, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this](double v){
    if (ui->slTStress && ui->slTStress->value()!=static_cast<int>(v))
      ui->slTStress->setValue(static_cast<int>(v));
    armStressDebounce_();
  });

  if (ui->cbStressAutoRun) {
    ui->cbStressAutoRun->setChecked(true);
    connect(ui->cbStressAutoRun, &QCheckBox::toggled, this, [this](bool on){
      if (on) armStressDebounce_(0);
    });
  }
  // === Vega Stress (%) ===
  // slider 0..100 → spin (double) and auto-run
  if (ui->slVegaStress) {
    connect(ui->slVegaStress, &QSlider::valueChanged, this, [this](int s){
      if (ui->sbVegaStress && std::lround(ui->sbVegaStress->value()) != s)
        ui->sbVegaStress->setValue(static_cast<double>(s));
      qDebug() << "[UI] sbVegaStress changed ->" << s << "%";
      armStressDebounce_();   // <- trigger auto stress
    });
  }

  // spin → slider and auto-run
  if (ui->sbVegaStress) {
    connect(ui->sbVegaStress, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, [this](double v){
      const int s = static_cast<int>(std::lround(v));
      if (ui->slVegaStress && ui->slVegaStress->value() != s)
        ui->slVegaStress->setValue(s);
      qDebug() << "[UI] slVegaStress changed ->" << v << "%";
      armStressDebounce_();   // <- trigger auto stress
    });
  }

  // If “Use calibrated surface” changes, re-run too
  if (ui->chkUseSmile) {
    connect(ui->chkUseSmile, &QCheckBox::toggled, this, [this](bool){
      armStressDebounce_();
    });
  }
}

void MainWindow::armStressDebounce_(int ms) {
  const bool on = (ui->cbStressAutoRun ? ui->cbStressAutoRun->isChecked() : true);
  if (!on) return;
  if (stressBusy_) { stressPending_ = true; return; }
  stressDebounce_->start(ms >= 0 ? ms : 0);
}

// ========================= Stress =========================
void MainWindow::setupStressChart() {
  if (stressChartView_) return;
  using namespace QtCharts;

  stressChart_ = new QChart();
  stressChart_->setTitle("Stress P&L – MC vs Taylor");
  stressChart_->setMargins(QMargins(8,8,18,12));
  stressChart_->legend()->setVisible(false);
  stressChart_->setAnimationOptions(QChart::SeriesAnimations);

  // Séries empilées Taylor (parent = chart / série)
  stressTaylorSeries_ = new QStackedBarSeries(stressChart_);
  auto* sDelta = new QBarSet("Δ·ΔS",     stressTaylorSeries_); *sDelta << 0.0 << 0.0;
  auto* sGamma = new QBarSet("½Γ(ΔS)^2", stressTaylorSeries_); *sGamma << 0.0 << 0.0;
  auto* sVega  = new QBarSet("Vega·Δσ",  stressTaylorSeries_); *sVega  << 0.0 << 0.0;
  auto* sRho   = new QBarSet("Rho·Δr",   stressTaylorSeries_); *sRho   << 0.0 << 0.0;

  // IMPORTANT : on les ajoute à la série (le parent QObject ne suffit pas)
  stressTaylorSeries_->append({sDelta, sGamma, sVega, sRho});

  sDelta->setColor(QColor("#42A5F5"));
  sGamma->setColor(QColor("#AB47BC"));
  sVega ->setColor(QColor("#FFB300"));
  sRho  ->setColor(QColor("#26A69A"));

  stressTaylorSeries_->setBarWidth(0.42);
  stressTaylorSeries_->setLabelsVisible(true);
  stressTaylorSeries_->setLabelsPosition(QAbstractBarSeries::LabelsOutsideEnd);
  stressTaylorSeries_->setLabelsFormat("@value");

  // Série MC (parent = chart / série)
  stressMCSeries_ = new QBarSeries(stressChart_);
  auto* sMC = new QBarSet("P&L (MC)", stressMCSeries_); *sMC << 0.0 << 0.0;
  sMC->setColor(QColor("#455A64"));
  stressMCSeries_->setBarWidth(0.42);
  stressMCSeries_->setLabelsVisible(true);
  stressMCSeries_->setLabelsPosition(QAbstractBarSeries::LabelsOutsideEnd);
  stressMCSeries_->setLabelsFormat("@value");

  // IMPORTANT : on ajoute le set MC à la série
  stressMCSeries_->append(sMC);

  stressChart_->addSeries(stressTaylorSeries_);
  stressChart_->addSeries(stressMCSeries_);

  // Axes (parent = chart)
  stressAxisX_ = new QBarCategoryAxis(stressChart_);
  stressAxisX_->append(QStringList() << "Taylor" << "MC");
  stressChart_->addAxis(stressAxisX_, Qt::AlignBottom);
  stressTaylorSeries_->attachAxis(stressAxisX_);
  stressMCSeries_->attachAxis(stressAxisX_);

  stressAxisY_ = new QValueAxis(stressChart_);
  stressAxisY_->setTitleText("P&L");
  stressAxisY_->setLabelFormat("%.2f");
  stressAxisY_->setMinorTickCount(1);
  stressAxisY_->setGridLineVisible(true);
  stressChart_->addAxis(stressAxisY_, Qt::AlignLeft);
  stressTaylorSeries_->attachAxis(stressAxisY_);
  stressMCSeries_->attachAxis(stressAxisY_);

  // Ligne zéro (parent = chart)
  auto* zeroLine = new QLineSeries(stressChart_);
  zeroLine->setUseOpenGL(false);
  zeroLine->setPen(QPen(QColor("#9E9E9E"), 1, Qt::DashLine));
  zeroLine->append(-0.4, 0.0);
  zeroLine->append( 1.4, 0.0);
  stressChart_->addSeries(zeroLine);
  zeroLine->attachAxis(stressAxisX_);
  zeroLine->attachAxis(stressAxisY_);

  // View parentée au container
  stressChartView_ = new QChartView(stressChart_, ui->stressChartContainer);
  stressChartView_->setRenderHint(QPainter::Antialiasing);
  stressChartView_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  stressChartView_->setMinimumHeight(280);

  auto* v = new QVBoxLayout(ui->stressChartContainer);
  v->setContentsMargins(0,0,0,0);
  v->addWidget(stressChartView_);
}

void MainWindow::resetStressChart() {
  if (!stressTaylorSeries_ || !stressMCSeries_) return;

  // Normalise chaque set à exactement 2 valeurs [Taylor, MC] puis remets à zéro
  auto normalizeAndZero = [](QtCharts::QBarSet* s){
    if (!s) return;
    const int n = s->count();
    if (n < 2) { s->remove(0, n); s->append(0.0); s->append(0.0); }
    else if (n > 2) { s->remove(2, n - 2); }
    // maintenant count()==2
    s->replace(0, 0.0);
    s->replace(1, 0.0);
  };

  for (auto* set : stressTaylorSeries_->barSets()) normalizeAndZero(set);
  for (auto* set : stressMCSeries_->barSets())     normalizeAndZero(set);

  if (stressAxisY_) {
    stressAxisY_->setRange(-1.0, 1.0);
    stressAxisY_->applyNiceNumbers();
  }
}

void MainWindow::updateStressChart(double pnlMC,
                                   double d_contrib, double g_contrib,
                                   double v_contrib, double r_contrib)
{
  if (!stressTaylorSeries_ || !stressMCSeries_) return;

  // Vérifie la présence des sets (évite tout out-of-range)
  if (stressTaylorSeries_->barSets().size() < 4 || stressMCSeries_->barSets().size() < 1) {
    qWarning() << "Stress chart not initialized: missing bar sets";
    return;
  }

  auto* sDelta = stressTaylorSeries_->barSets().at(0);
  auto* sGamma = stressTaylorSeries_->barSets().at(1);
  auto* sVega  = stressTaylorSeries_->barSets().at(2);
  auto* sRho   = stressTaylorSeries_->barSets().at(3);
  auto* sMC    = stressMCSeries_->barSets().at(0);

  // Normalise chaque set à exactement 2 valeurs [Taylor, MC]
  auto ensure2 = [](QtCharts::QBarSet* s){
    if (!s) return;
    const int n = s->count();
    if (n < 2) { s->remove(0, n); s->append(0.0); s->append(0.0); }
    else if (n > 2) { s->remove(2, n - 2); }
  };
  ensure2(sDelta); ensure2(sGamma); ensure2(sVega); ensure2(sRho); ensure2(sMC);

  // Sanity numérique
  auto sane = [](double x){ return std::isfinite(x) ? x : 0.0; };
  pnlMC     = sane(pnlMC);
  d_contrib = sane(d_contrib);
  g_contrib = sane(g_contrib);
  v_contrib = sane(v_contrib);
  r_contrib = sane(r_contrib);

  // Remplacements atomiques (pas de remove/append qui crée des états transitoires)
  sDelta->replace(0, d_contrib); sDelta->replace(1, 0.0);
  sGamma->replace(0, g_contrib); sGamma->replace(1, 0.0);
  sVega ->replace(0, v_contrib); sVega ->replace(1, 0.0);
  sRho  ->replace(0, r_contrib); sRho  ->replace(1, 0.0);

  sMC->replace(0, 0.0);
  sMC->replace(1, pnlMC);

  const double taylor = d_contrib + g_contrib + v_contrib + r_contrib;
  double lo = std::min({taylor, pnlMC, 0.0});
  double hi = std::max({taylor, pnlMC, 0.0});
  if (!std::isfinite(lo) || !std::isfinite(hi)) { lo = -1.0; hi = 1.0; }

  const double pad = std::max(1.0, 0.08 * std::max(std::abs(lo), std::abs(hi)));
  if (stressAxisY_) {
    stressAxisY_->setRange(lo - pad, hi + pad);
    stressAxisY_->applyNiceNumbers();  // arrondit les graduations
  }

  // (optionnel) Si l’axe X a été modifié ailleurs, on rétablit 2 catégories
  if (stressAxisX_ && stressAxisX_->categories().size() != 2) {
    stressAxisX_->clear();
    stressAxisX_->append(QStringList() << "Taylor" << "MC");
  }

  // Si tu décides d’afficher la légende native :
  // stressChart_->legend()->setVisible(true);
  stressChart_->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
  stressChart_->legend()->setAlignment(Qt::AlignBottom);
}

void MainWindow::setupStressLegend() {
  if (!ui->legendLayoutStress) return;

  // nettoie
  QLayoutItem* it;
  while ((it = ui->legendLayoutStress->takeAt(0))) {
    delete it->widget(); delete it;
  }

  auto add = [&](const QString& label, const QColor& c){
    auto* sw = new QLabel; sw->setFixedSize(18,12);
    QPixmap pm(18,12); pm.fill(c); sw->setPixmap(pm);
    auto* row = new QHBoxLayout; row->setSpacing(6);
    row->addWidget(sw); row->addWidget(new QLabel(label)); row->addStretch();
    auto* w = new QWidget; w->setLayout(row);
    ui->legendLayoutStress->addWidget(w);
  };

  add("Δ·ΔS",      QColor("#42A5F5"));
  add("½Γ(ΔS)^2",  QColor("#AB47BC"));
  add("Vega·Δσ",   QColor("#FFB300"));
  add("Rho·Δr",    QColor("#26A69A"));
  add("P&L (MC)",  QColor("#455A64"));
  ui->legendLayoutStress->addStretch();
}


/* 
===================================================================
============================== Test ===============================
===================================================================
*/

void MainWindow::testsClearTable() const {
  if (!ui->twTests) return;
  ui->twTests->clear();
  ui->twTests->setColumnCount(5);
  ui->twTests->setHorizontalHeaderLabels({"Test","Metric","Value","Threshold","Result"});
  ui->twTests->setRowCount(0);
}

void MainWindow::testsLog(const QString& line) const {
  if (ui->pteLog) ui->pteLog->appendPlainText(line);
}

double MainWindow::slopeLeastSquares(const std::vector<double>& x, const std::vector<double>& y) {
  const int n = (int)std::min(x.size(), y.size());
  if (n < 2) return std::numeric_limits<double>::quiet_NaN();
  double sx=0, sy=0, sxx=0, sxy=0;
  for (int i=0;i<n;++i) { sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; }
  const double d = n*sxx - sx*sx;
  if (std::fabs(d) < 1e-18) return std::numeric_limits<double>::quiet_NaN();
  return (n*sxy - sx*sy) / d; // pente
}

void MainWindow::runAllTests() {
  testsClearTable();
  testsLog("=== Validation start ===");
  if (mcWorker_) QObject::disconnect(mcWorker_, nullptr, this, nullptr); // sécurité
  if (ui->cbCRN) ui->cbCRN->setChecked(true);
  runTestPricingBS();
}

void MainWindow::runTestPricingBS() {
  testsLog("[Pricing vs BS]");

  const double S0=ui->sbS0->value(), r=ui->sbR->value(), q=ui->sbQ->value();
  const double sigma=ui->sbSigma->value(), K=ui->sbK->value(), T=ui->sbT->value();
  const bool isPut=(ui->cbType->currentText().toLower()=="put");

  rw::market::MarketData mkt(S0,r,q);
  rw::market::Instrument inst(K,T, isPut?rw::market::OptionType::Put:rw::market::OptionType::Call);

  const std::size_t nTarget=(std::size_t)ui->sbNTarget->value();
  const std::size_t batch  =(std::size_t)ui->sbBatch->value();
  const double      tol    =ui->dsbTol->value();
  const std::size_t nSteps =(std::size_t)ui->sbNSteps->value();
  const std::uint64_t seed =(std::uint64_t)ui->sbSeed->value();
  rw::config::McConfig cfg(nTarget,batch,tol,nSteps,seed);

  const double bs = bs_price(mkt, inst, sigma);

  startMcWorker();
  applySmileSettingsToWorker_();
  connect(mcWorker_, &gui::McWorker::finished, this,
          [=](double price,double se,double lo,double hi,std::size_t nEff,long long ms){
            const bool ok = (bs>=lo && bs<=hi);
            TestRow row;
            row.test="Pricing vs BS";
            row.metric="BS in CI95%";
            row.value=QString::asprintf("BS=%.6f CI=[%.6f, %.6f]", bs, lo, hi);
            row.threshold="inside";
            row.pass=ok;
            testsAppendRow(row);
            testsLog(QString::asprintf("  price=%.6f, se=%.3g, BS=%.6f -> %s",
                           price, se, bs, ok ? "PASS" : "FAIL"));

            stopMcWorker();
            // Enchaîne
            runTestConvergence();
          },
          Qt::UniqueConnection);

  QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, mkt, inst, cfg, sigma]{ w->runPricing(mkt, inst, cfg, sigma); },
      Qt::QueuedConnection);
}

void MainWindow::runTestConvergence() {
  testsLog("[Convergence 1/sqrt(N)]");
  convN_.clear(); convHalf_.clear();

  const double S0=ui->sbS0->value(), r=ui->sbR->value(), q=ui->sbQ->value();
  const double sigma=ui->sbSigma->value(), K=ui->sbK->value(), T=ui->sbT->value();
  const bool isPut=(ui->cbType->currentText().toLower()=="put");

  rw::market::MarketData mkt(S0,r,q);
  rw::market::Instrument inst(K,T, isPut?rw::market::OptionType::Put:rw::market::OptionType::Call);

  const std::size_t nTarget=(std::size_t)ui->sbNTarget->value();
  const std::size_t batch  =(std::size_t)ui->sbBatch->value();
  const std::size_t nSteps =(std::size_t)ui->sbNSteps->value();
  const std::uint64_t seed =(std::uint64_t)ui->sbSeed->value();
  rw::config::McConfig cfg(nTarget,batch,-1.0,nSteps,seed);

  startMcWorker();
  applySmileSettingsToWorker_();

  connect(mcWorker_, &gui::McWorker::progress, this,
          [this](std::size_t N, double /*mean*/, double half){
            if (N>0 && half>1e-16) { convN_.push_back(std::log((double)N));
                                     convHalf_.push_back(std::log(half)); }
          }, Qt::UniqueConnection);

  connect(mcWorker_, &gui::McWorker::finished, this,
          [this](double /*p*/, double /*se*/, double /*lo*/, double /*hi*/,
                 std::size_t, long long){
            const double b = slopeLeastSquares(convN_, convHalf_); // attendu ~ -0.5
            const bool ok = std::isfinite(b) && (b>=-0.60 && b<=-0.40);
            TestRow row;
            row.test="Convergence";
            row.metric="slope(log half) vs log N";
            row.value=QString::asprintf("%.3f", b);
            row.threshold="[-0.60,-0.40]";
            row.pass=ok;
            testsAppendRow(row);
            testsLog(QString::asprintf("  slope=%.3f -> %s", b, ok ? "PASS" : "FAIL"));
            stopMcWorker();
            // Enchaîne
            runTestGreeksBS();
          },
          Qt::UniqueConnection);

  QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, mkt, inst, cfg, sigma]{ w->runPricing(mkt, inst, cfg, sigma); },
      Qt::QueuedConnection);
}

void MainWindow::runTestGreeksBS() {
  testsLog("Test: Greeks vs Black–Scholes (z-score ≤ 2)");

  // 1) État UI
  const double S0    = ui->sbS0->value();
  const double r     = ui->sbR->value();
  const double q     = ui->sbQ->value();
  const double sigma = ui->sbSigma->value();
  const double K     = ui->sbK->value();
  const double T     = ui->sbT->value();
  const bool   isPut = (ui->cbType->currentText().toLower() == "put");

  rw::market::MarketData mkt(S0, r, q);
  rw::market::OptionType typ = isPut ? rw::market::OptionType::Put
                                     : rw::market::OptionType::Call;
  rw::market::Instrument inst(K, T, typ);

  // 2) MC/config
  const std::size_t nTarget = static_cast<std::size_t>(ui->sbNTarget->value());
  const std::size_t batch   = static_cast<std::size_t>(ui->sbBatch->value());
  const double      tol     = ui->dsbTol->value();
  const std::size_t nSteps  = static_cast<std::size_t>(ui->sbNSteps->value());
  const std::uint64_t seed  = static_cast<std::uint64_t>(ui->sbSeed->value());
  rw::config::McConfig mc(nTarget, batch, tol, nSteps, seed);

  rw::config::GreeksConfig gk(rw::config::Greek::Delta);
  gk.bump_rel_S0    = ui->dsbBumpRelS0->value();
  gk.bump_abs_sigma = ui->dsbBumpAbsSigma->value();
  gk.bump_abs_r     = ui->dsbBumpAbsR->value();
  gk.bump_abs_T     = ui->dsbBumpAbsT->value();
  gk.use_antithetic = ui->cbAnti->isChecked();
  gk.use_crn        = true;            // stabilise la campagne
  gk.n_steps        = nSteps;

  // ⚠️ clamp du bump en sigma : ≤ 5% de σ, et au plus 1e-3
  {
    const double cap = std::max(1e-4, 0.05 * std::max(1e-8, sigma)); // 5% de σ, min 1e-4
    if (gk.bump_abs_sigma > cap) {
      testsLog(QString::asprintf("  clamp bump σ: %.6f -> %.6f", gk.bump_abs_sigma, cap));
      gk.bump_abs_sigma = cap;
    }
  }

  // Pour le test : on force la méthode à "lrm" (Vega via LRM, Γ via Δ PW)
  QString method = "lrm";

  // 3) État de campagne
  struct State {
    int reps{10};
    int done{0};
    std::vector<double> D,V,G,R,T;
    QMetaObject::Connection conn;
  };
  auto st = std::make_shared<State>();
  st->D.reserve(st->reps);
  st->V.reserve(st->reps);
  st->G.reserve(st->reps);
  st->R.reserve(st->reps);
  st->T.reserve(st->reps);

  // 4) Lancement (seed décalé à chaque rép)
  auto runOnce = [this, mkt, inst, mc, gk, sigma, method, st]() {
    rw::config::McConfig mc2 = mc;
    mc2.seed = mc.seed + static_cast<std::uint64_t>(st->done);
    QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, mkt, inst, mc2, gk, sigma, method]() {
        w->runGreeks(mkt, inst, mc2, gk, sigma, method);
      },
      Qt::QueuedConnection
    );
  };

  startMcWorker();
  applySmileSettingsToWorker_();
  st->conn = connect(mcWorker_, &gui::McWorker::greeksFinished, this,
    [this, st, mkt, inst, sigma, runOnce](double d,double v,double g,double r,double t,long long){
      st->D.push_back(d); st->V.push_back(v); st->G.push_back(g); st->R.push_back(r); st->T.push_back(t);
      st->done++;
      if (st->done < st->reps) { runOnce(); return; }

      disconnect(st->conn);

      auto mean = [](const std::vector<double>& x){
        if (x.empty()) return 0.0;
        return std::accumulate(x.begin(), x.end(), 0.0) / double(x.size());
      };
      double Dm = mean(st->D), Vm = mean(st->V), Gm = mean(st->G),
             Rm = mean(st->R), Tm = mean(st->T);

      auto se = [this](const std::vector<double>& v) {
        const std::size_t n = v.size();
        if (n < 2) return 0.0;
        return std::sqrt(sampleVariance(v, /*unbiased*/true) / double(n));
      };
      double Dse = se(st->D), Vse = se(st->V), Gse = se(st->G),
             Rse = se(st->R), Tse = se(st->T);

      // 5) Réf BS
      const auto ga = bs_greeks(mkt, inst, sigma);
      const double eps = 1e-12;

      // 6) Auto-normalisation d’unités
      auto snap_and_scale = [](double mc, double bs, double candidate, double tol_rel=0.35){
        if (std::fabs(mc) < 1e-16 || std::fabs(bs) < 1e-16) return mc;
        const double k = std::fabs(bs / mc);
        if (std::fabs(k - candidate) <= tol_rel * candidate) return mc * candidate;
        return mc;
      };
      Vm = snap_and_scale(Vm, ga.vega,   100.0);   // vega par %
      Rm = snap_and_scale(Rm, ga.rho,  10000.0);   // rho par bp
      Tm = snap_and_scale(Tm, ga.theta, 365.0);    // theta par jour

      auto safe_ratio = [eps](double mc, double bs){ return (std::fabs(bs) < eps ? 0.0 : mc/bs); };
      testsLog(QString::asprintf(
        "  Ratios MC/BS après normalisation: Δ=%.3f  V=%.3f  Γ=%.3f  ρ=%.3f  Θ=%.3f",
        safe_ratio(Dm,ga.delta), safe_ratio(Vm,ga.vega), safe_ratio(Gm,ga.gamma),
        safe_ratio(Rm,ga.rho),   safe_ratio(Tm,ga.theta)));

      // 7) z-scores
      auto z = [](double m, double ref, double se){ return (se <= 0.0) ? 0.0 : ((m - ref) / se); };
      const double zD = z(Dm, ga.delta, Dse);
      const double zV = z(Vm, ga.vega,  Vse);
      const double zG = z(Gm, ga.gamma, Gse);
      const double zR = z(Rm, ga.rho,   Rse);
      const double zT = z(Tm, ga.theta, Tse);

      const double zMax = std::max({std::fabs(zD), std::fabs(zV), std::fabs(zG),
                                    std::fabs(zR), std::fabs(zT)});
      const bool pass = (zMax <= 2.0);

      testsAppendRow({ "Greeks vs BS", "max |z|",
                       QString::asprintf("%.3f", zMax), "≤ 2.0",
                       pass ? "PASS" : "FAIL", pass });

      testsLog(QString::asprintf(
        "  MC: Δ=%.6f V=%.6f Γ=%.6f ρ=%.6f Θ=%.6f | BS: Δ=%.6f V=%.6f Γ=%.6f ρ=%.6f Θ=%.6f",
        Dm,Vm,Gm,Rm,Tm, ga.delta,ga.vega,ga.gamma,ga.rho,ga.theta));
      testsLog(QString::asprintf(
        "  z:  Δ=%.2f  V=%.2f  Γ=%.2f  ρ=%.2f  Θ=%.2f  -> %s",
        zD,zV,zG,zR,zT, pass?"PASS":"FAIL"));

      stopMcWorker();
      runTestStressCRN();
    },
    Qt::UniqueConnection);

  runOnce();
}

void MainWindow::runTestStressCRN() {
  testsLog("[Stress CRN: Var(P&L) << Var(price)]");
  if (!baselineSet_) onSnapBaseline();

  // Chocs UI -> unités internes
  const double dS_rel   = (ui->sbS0Stress    ? ui->sbS0Stress->value()    : 0.0) / 100.0;
  const double dSigma   =  ui->sbSigmaStress ? ui->sbSigmaStress->value() : 0.0;
  const double dR_abs   = (ui->sbRStress     ? ui->sbRStress->value()     : 0.0) / 10000.0;
  const double dT_years = (ui->sbTStress     ? ui->sbTStress->value()     : 0.0) / 365.0;

  // Si l'utilisateur n'a mis aucun choc, applique un stress "canonique"
  const double stress_norm = std::fabs(dS_rel) + std::fabs(dSigma)
                          + std::fabs(dR_abs) + std::fabs(dT_years);
  if (stress_norm < 1e-12) {
    // Defaults raisonnables pour un test non trivial
    double dS_rel   = 0.05;           // +5%
    double dSigma   = 0.01;           // +1 pt de vol
    double dR_abs   = -0.0025;        // -25 bps
    double dT_years = 30.0 / 365.0;   // +30 jours
    testsLog("  (no UI stress) using defaults: dS=+5%, dσ=+0.01, dr=-25bps, dT=+30d");
  }

  // Paramètres MC (on ne stocke PAS McConfig, seulement les scalaires)
  const std::size_t nTarget = static_cast<std::size_t>(ui->sbNTarget->value());
  const std::size_t batch   = static_cast<std::size_t>(ui->sbBatch->value());
  const std::size_t nSteps  = static_cast<std::size_t>(ui->sbNSteps->value());
  const std::uint64_t seed0 = static_cast<std::uint64_t>(ui->sbSeed->value());

  struct State {
    int Krep = 20;
    int rep  = 0;
    std::vector<double> pnl, px;
    double dS_rel{}, dSigma{}, dR_abs{}, dT_years{};
    std::size_t nTarget{}, batch{}, nSteps{};
    std::uint64_t seed0{};
  };
  auto st = std::make_shared<State>();
  st->dS_rel   = dS_rel;
  st->dSigma   = dSigma;
  st->dR_abs   = dR_abs;
  st->dT_years = dT_years;
  st->nTarget  = nTarget;
  st->batch    = batch;
  st->nSteps   = nSteps;
  st->seed0    = seed0;
  st->pnl.reserve(st->Krep);
  st->px .reserve(st->Krep);

  startMcWorker();
  applySmileSettingsToWorker_();

  auto launch = [this, st]() {
    rw::config::McConfig cfgR(
      st->nTarget,
      st->batch,
      -1.0,               // tolérance off
      st->nSteps,
      st->seed0 + 2000u + static_cast<std::uint64_t>(st->rep)
    );
    QMetaObject::invokeMethod(
      mcWorker_,
      [w=mcWorker_, this, st, cfgR]() {
        w->runStress(*baseMkt_, *baseInst_, cfgR, baseSigma_,
                     st->dS_rel, st->dSigma, st->dR_abs, st->dT_years);
      },
      Qt::QueuedConnection
    );
  };

  connect(mcWorker_, &gui::McWorker::stressFinished, this,
          [this, st, launch](double /*baseP*/, double stressP, double pnlMC, long long){
            st->pnl.push_back(pnlMC);
            st->px .push_back(stressP);

            if (++st->rep < st->Krep) { launch(); return; }

            if (st->pnl.size() < 2 || st->px.size() < 2) {
              TestRow row; row.test="Stress (CRN)"; row.metric="samples";
              row.value=QString::number(static_cast<int>(std::min(st->pnl.size(), st->px.size())));
              row.threshold="≥ 2"; row.pass=false; row.result="FAIL";
              testsAppendRow(row);
              testsLog("  not enough samples -> FAIL");
              stopMcWorker();
              runTestSaveReload();
              return;
            }

            const double varPL    = sampleVariance(st->pnl, /*unbiased*/true);
            const double varPrice = sampleVariance(st->px , /*unbiased*/true);

            double ratio = std::numeric_limits<double>::infinity();
            if (varPrice > 0.0) ratio = varPL / varPrice;

            const bool ok = std::isfinite(ratio) && (ratio <= 0.20);

            TestRow row;
            row.test      = "Stress (CRN)";
            row.metric    = "Var(P&L)/Var(price)";
            row.value     = QString::asprintf("%.3e (VarPL=%.3e, VarPrice=%.3e)", ratio, varPL, varPrice);
            row.threshold = "≤ 0.20";
            row.pass      = ok;
            row.result    = ok ? "PASS" : "FAIL";
            testsAppendRow(row);
            testsLog(QString::asprintf("  ratio=%.3e -> %s", ratio, ok ? "PASS" : "FAIL"));

            // petit résumé des échantillons
            auto summarize = [](const std::vector<double>& v){
              if (v.empty()) return QString("n=0");
              auto mm = std::minmax_element(v.begin(), v.end());
              return QString::asprintf("n=%zu, min=%.6g, max=%.6g",
                                       v.size(), *mm.first, *mm.second);
            };
            testsLog("  P&L sample:   " + summarize(st->pnl));
            testsLog("  Price sample: " + summarize(st->px));

            stopMcWorker();
            runTestSaveReload();
          },
          Qt::UniqueConnection);

  launch();
}

void MainWindow::runTestSaveReload() {
  testsLog("[Save / Load JSON]");

  // Référence de pricing (à partir de l’UI)
  const double S0=ui->sbS0->value(), r=ui->sbR->value(), q=ui->sbQ->value();
  const double sigma=ui->sbSigma->value(), K=ui->sbK->value(), T=ui->sbT->value();
  const bool isPut=(ui->cbType->currentText().toLower()=="put");

  rw::market::MarketData mkt(S0,r,q);
  rw::market::Instrument inst(K,T, isPut?rw::market::OptionType::Put:rw::market::OptionType::Call);

  const std::size_t nTarget=(std::size_t)ui->sbNTarget->value();
  const std::size_t batch  =(std::size_t)ui->sbBatch->value();
  const double      tol    =ui->dsbTol->value();
  const std::size_t nSteps =(std::size_t)ui->sbNSteps->value();
  const std::uint64_t seed =(std::uint64_t)ui->sbSeed->value();
  rw::config::McConfig cfg(nTarget,batch,tol,nSteps,seed);

  // Snapshot JSON
  QJsonObject snapshot = makeProjectJson();

  // 1) Prix 1
  startMcWorker();
  applySmileSettingsToWorker_();
  auto conn1 = std::make_shared<QMetaObject::Connection>();
  *conn1 = connect(mcWorker_, &gui::McWorker::finished, this,
    [=](double p0,double,double,double,std::size_t,long long){
      QObject::disconnect(*conn1);
      stopMcWorker();                // stop premier worker

      // 2) Reload JSON
      loadProjectJson(snapshot);

      // 3) Prix 2 avec NOUVEAU worker
      startMcWorker();
      applySmileSettingsToWorker_();
      auto conn2 = std::make_shared<QMetaObject::Connection>();
      *conn2 = connect(mcWorker_, &gui::McWorker::finished, this,
        [=](double p1,double,double,double,std::size_t,long long){
          QObject::disconnect(*conn2);
          const double diff = std::fabs(p1 - p0);
          const bool ok = (diff <= 1e-9);

          testsAppendRow({ "Save/Load", "|price2 - price1|",
                           QString::asprintf("%.3g", diff), "≤ 1e-9",
                           ok ? "PASS" : "FAIL", ok });
          testsLog(QString::asprintf("  Δ=%.3g -> %s", diff, ok?"PASS":"FAIL"));

          stopMcWorker();
          testsLog("=== Validation end ===");
        },
        Qt::UniqueConnection);

      QMetaObject::invokeMethod(
        mcWorker_,
        [w=mcWorker_, mkt, inst, cfg, sigma]{ w->runPricing(mkt, inst, cfg, sigma); },
        Qt::QueuedConnection);
    },
    Qt::UniqueConnection);

  QMetaObject::invokeMethod(
    mcWorker_,
    [w=mcWorker_, mkt, inst, cfg, sigma]{ w->runPricing(mkt, inst, cfg, sigma); },
    Qt::QueuedConnection);
}


// ===== Tableau & log =====

void MainWindow::testsAppendRow(const TestRow& r) const {
  if (!ui->twTests) return;
  int row = ui->twTests->rowCount();
  ui->twTests->insertRow(row);
  auto put = [&](int c, const QString& s){
    auto* it = new QTableWidgetItem(s);
    if (c==4) {
      it->setForeground(r.pass ? QBrush(QColor("#1B5E20")) : QBrush(QColor("#B71C1C")));
      it->setText(r.pass ? "PASS" : "FAIL");
      it->setTextAlignment(Qt::AlignCenter);
    }
    ui->twTests->setItem(row, c, it);
  };
  put(0,r.test); put(1,r.metric); put(2,r.value); put(3,r.threshold); put(4,r.result);
  ui->twTests->resizeColumnsToContents();
}


// ===== Stats =====

double MainWindow::sampleVariance(const std::vector<double>& x, bool unbiased) {
  const std::size_t n = x.size();
  if (n < 2) return 0.0;  // évite toute lecture hors borne

  // Welford: stable numériquement, pas d'accès x[0] si vide
  double mean = 0.0;
  double M2   = 0.0;
  std::size_t k = 0;
  for (double v : x) {
    ++k;
    const double delta  = v - mean;
    mean += delta / static_cast<double>(k);
    const double delta2 = v - mean;
    M2 += delta * delta2;
  }
  return unbiased ? (M2 / static_cast<double>(n - 1))
                  : (M2 / static_cast<double>(n));
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        Calibration          %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

void MainWindow::wireCalibration_() {
  // Thread + worker
  calibThread_ = new QThread(this);
  calib_ = new gui::CalibWorker();
  calib_->moveToThread(calibThread_);
  connect(calibThread_, &QThread::finished, calib_, &QObject::deleteLater);
  calibThread_->start();

  // UI -> worker (queued)
  connect(ui->btnCalibOpen,   &QPushButton::clicked, this, &MainWindow::onCalibOpenCsv);
  connect(ui->btnCalibInvert,    &QPushButton::clicked, this, &MainWindow::onCalibInvertIv);
  connect(ui->btnCalibSmile,&QPushButton::clicked, this, &MainWindow::onCalibBuildSmile);
  connect(ui->btnCalibFitGlobal, &QPushButton::clicked, this, &MainWindow::onCalibFitGlobal);
  connect(ui->btnCalibRunQc,     &QPushButton::clicked, this, &MainWindow::onCalibRunQc);
  connect(ui->btnCalibExport,    &QPushButton::clicked, this, &MainWindow::onCalibExportCsv);
  connect(ui->btnCalibSaveProject,&QPushButton::clicked,this, &MainWindow::onCalibSaveProject);
  connect(ui->cbCalibT, qOverload<int>(&QComboBox::currentIndexChanged),
          this, &MainWindow::onCalibSelectT);

  // worker -> UI
  connect(calib_, &gui::CalibWorker::message,       this, &MainWindow::onCalibMessage);
  connect(calib_, &gui::CalibWorker::progress,      this, &MainWindow::onCalibProgress);
  connect(calib_, &gui::CalibWorker::loadedMarket,  this, &MainWindow::onCalibLoaded);
  connect(calib_, &gui::CalibWorker::inverted,      this, &MainWindow::onCalibInverted);
  connect(calib_, &gui::CalibWorker::smileBuilt,    this, &MainWindow::onCalibSmileBuilt);
  connect(calib_, &gui::CalibWorker::globalFitted,  this, &MainWindow::onCalibGlobalFitted);
  connect(calib_, &gui::CalibWorker::qcDone,        this, &MainWindow::onCalibQcDone);
  connect(calib_, &gui::CalibWorker::exportedCsv,   this, &MainWindow::onCalibExported);
  connect(calib_, &gui::CalibWorker::projectSaved,  this, &MainWindow::onCalibProjectSaved);
  connect(calib_, &gui::CalibWorker::failed,        this, &MainWindow::onCalibFailed);
}

void MainWindow::onCalibOpenCsv() {
  const QString path = QFileDialog::getOpenFileName(
      this, "Choisir un CSV marché", "data/market_samples", "CSV (*.csv)");
  if (path.isEmpty()) return;
  QMetaObject::invokeMethod(calib_, "loadCsv",
                            Qt::QueuedConnection,
                            Q_ARG(QString, path));
}

void MainWindow::onCalibInvertIv() {
  qDebug() << "[UI] request invertIv()";
  QMetaObject::invokeMethod(calib_, "invertIv", Qt::QueuedConnection);
}

void MainWindow::onCalibSmileBuilt(const QVector<double>& maturities)
{
  qDebug() << "[UI] Smile built: slices=" << maturities.size() << " mats=" << maturities;
  smileCalibCache_ = calib_->smile();

  ui->cbCalibT->clear();
  for (double T : maturities) ui->cbCalibT->addItem(QString::number(T));
  if (ui->cbCalibT->count()>0) ui->cbCalibT->setCurrentIndex(0);

  repaintSmile_();
  repaintHeatmap_();
  finishCalibStage_(true, "Smile construit");
}

void MainWindow::onCalibFitGlobal() {
  QMetaObject::invokeMethod(calib_, "fitGlobal",
                            Qt::QueuedConnection,
                            Q_ARG(double, 0.20));
}

void MainWindow::onCalibRunQc() {
  // tick=0, tol=1.0 pour "only-consistent", pas de price-from-iv
  QMetaObject::invokeMethod(calib_, "runQc",
                            Qt::QueuedConnection,
                            Q_ARG(double, 0.0),
                            Q_ARG(double, -1.0),
                            Q_ARG(bool,   true));
}

void MainWindow::onCalibExportCsv() {
  QMetaObject::invokeMethod(calib_, "exportSurfaceCsv",
                            Qt::QueuedConnection,
                            Q_ARG(QString, QString("data/iv_exports")));
}

void MainWindow::onCalibSaveProject() {
  QMetaObject::invokeMethod(calib_, "saveProjectJson", Qt::QueuedConnection);
}


void MainWindow::onCalibMessage(const QString& m) {
  statusBar()->showMessage(m, 3000);
}

void MainWindow::onCalibProgress(const QString& stage, int cur, int tot){
  ui->lblCalibMeta->setText(stage);
  if (cur == 0) startCalibStage_(stage, tot);       // (ré)initialise + rouge
  else          updateCalibStage_(cur, tot);        // met à jour + ETA

  if (tot>0) { ui->calibProgress->setRange(0,tot); ui->calibProgress->setValue(cur); }
  else { ui->calibProgress->setRange(0,0); }
}

void MainWindow::onCalibLoaded(size_t nRows, double S0, double r, double q, const QString& und)
{
  msCalibCache_ = calib_->market();
  ui->lblCalibMeta->setText(
      QString("underlying=%1  S0=%2  r=%3  q=%4  rows=%5")
      .arg(und).arg(S0).arg(r).arg(q).arg((qulonglong)nRows));

  // preview table
  auto& rows = msCalibCache_.rows;
  ui->tblCalibPreview->setRowCount((int)rows.size());
  for (int i=0; i<(int)rows.size(); ++i) {
    const auto& rr = rows[i];
    auto cell = [&](int c, const QString& s){ auto* it = new QTableWidgetItem(s); it->setFlags(it->flags() ^ Qt::ItemIsEditable); ui->tblCalibPreview->setItem(i,c,it); };
    cell(0, QString::number(rr.K));
    cell(1, QString::number(rr.T));
    cell(2, rr.is_call ? "C" : "P");
    cell(3, QString::number(rr.mid_price));
    cell(4, QString::number(rr.iv_mid));
    cell(5, QString::number(rr.bid));
    cell(6, QString::number(rr.ask));
    cell(7, ""); // note
  }
  ui->calibProgress->setRange(0,1); ui->calibProgress->setValue(1);
  finishCalibStage_(true, "CSV chargé");
}

void MainWindow::onCalibInverted(size_t /*priced*/, size_t /*ok*/, size_t /*bad*/)
{
  msCalibCache_ = calib_->market();
  // qDebug() << "[UI] IV inversion done: priced=" << priced << " ok=" << ok << " bad=" << bad;
  for (int i=0;i<(int)msCalibCache_.rows.size();++i)
    ui->tblCalibPreview->item(i,4)->setText(QString::number(msCalibCache_.rows[i].iv_mid));
  ui->calibProgress->setRange(0,1); ui->calibProgress->setValue(1);
  finishCalibStage_(true, "IV inversées");
}


void MainWindow::onCalibGlobalFitted(const rw::calib::CalibReport& rep, double sigma)
{
  statusBar()->showMessage(
      QString("σ*=%1  RMSE_price=%2  n=%3  %4")
      .arg(sigma).arg(rep.rmse_price).arg((qulonglong)rep.n)
      .arg(rep.converged ? "converged" : "not converged"), 4000);
    finishCalibStage_(true, "Calibration σ globale");
}

void MainWindow::onCalibQcDone(double rmse_price,
                               double rmse_iv,
                               size_t n,
                               const std::vector<rw::qc::ErrorRow>& rows){
  // 1) Mettre à jour le cache
  residCalibCache_ = rows;

  // 2) Logger un petit résumé
  const double first = rows.empty() ? 0.0 : rows.front().err_price;
  qDebug() << "[UI] QC reçu:"
           << "n=" << n
           << "rmse_price=" << rmse_price
           << "rmse_iv="    << rmse_iv
           << "first_err="  << first;

  // 3) Status bar + progression (vert = terminé)
  statusBar()->showMessage(
      QString("QC: n=%1  RMSE_price=%2  RMSE_iv=%3")
          .arg((qulonglong)n).arg(rmse_price).arg(rmse_iv),
      4000);
  finishCalibStage_(true, "QC terminé");   // ta fonction qui passe la barre en vert

  // 4) Repeindre le graphe de résidus
  repaintResiduals_();
  finishCalibStage_(true, "QC terminé");
}

void MainWindow::onCalibSelectT(int) { repaintSmile_(); }

void MainWindow::onCalibExported(const QString& path) {
  showFileToast_("Export CSV", path);
}

void MainWindow::onCalibProjectSaved(const QString& path) {
  showFileToast_("Projet sauvegardé", path);
}


void MainWindow::onCalibFailed(const QString& why) {
  statusBar()->showMessage(why, 5000);
  if (auto* p = findChild<QProgressBar*>("progCalib")) { p->setRange(0,1); p->setValue(0); }
  finishCalibStage_(false, "Erreur");
}

void MainWindow::repaintSmile_()
{
  if (ui->cbCalibT->count()==0) return;
  const int idx = ui->cbCalibT->currentIndex();
  if (idx < 0) return;
  const double T = ui->cbCalibT->itemText(idx).toDouble();
  updateSmileChartForT_(T);
}

void MainWindow::repaintHeatmap_()
{
  auto* lbl = ui->lblCalibHeatmap;             // ton QLabel dans le .ui
  if (!lbl) return;

  lbl->clear();
  if (smileCalibCache_.slices.empty()) return;
  if (msCalibCache_.S0 <= 0) return;

  // bornes K/T & IV
  double Tmin=1e300, Tmax=-1e300, kmin=1e300, kmax=-1e300, ivmin=1e300, ivmax=-1e300;
  for (const auto& s : smileCalibCache_.slices) {
    Tmin = std::min(Tmin, s.T);
    Tmax = std::max(Tmax, s.T);
    for (double x : s.x) {
      double K = msCalibCache_.S0 * std::exp(x);   // x = ln(K/S0)
      kmin = std::min(kmin, K);
      kmax = std::max(kmax, K);
    }
    for (double iv : s.y) {
      if (std::isfinite(iv)) { ivmin = std::min(ivmin, iv); ivmax = std::max(ivmax, iv); }
    }
  }
  if (!(Tmin < Tmax && kmin < kmax && ivmin < ivmax)) return;

  // image
  const int W = std::max(180, lbl->width());
  const int H = std::max(180, lbl->height());
  QImage img(W, H, QImage::Format_RGB32);
  img.fill(Qt::white);

  for (int y = 0; y < H; ++y) {
    const double Ty = Tmin + (1.0 - double(y) / (H - 1)) * (Tmax - Tmin);
    for (int x = 0; x < W; ++x) {
      const double Kx = std::exp(std::log(kmin) + (double(x) / (W - 1)) * (std::log(kmax) - std::log(kmin)));
      double iv = smileCalibCache_.iv_at(Ty, Kx);
      if (!std::isfinite(iv)) { img.setPixel(x, y, qRgb(240, 240, 240)); continue; }
      double t = (iv - ivmin) / (ivmax - ivmin);
      t = std::clamp(t, 0.0, 1.0);
      // simple colormap: blanc↔rouge
      int g = int(std::lround(255.0 * (1.0 - t)));
      img.setPixel(x, y, qRgb(255, g, g));
    }
  }
  qDebug() << "Heatmap range: ivmin=" << ivmin << " ivmax=" << ivmax
         << " Tmin=" << Tmin << " Tmax=" << Tmax
         << " kmin=" << kmin << " kmax=" << kmax;

  lbl->setPixmap(QPixmap::fromImage(img).scaled(lbl->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
}

void MainWindow::repaintResiduals_()
{
  std::vector<double> errs; errs.reserve(residCalibCache_.size());
  for (const auto& e : residCalibCache_) errs.push_back(e.err_price);
  updateResidChart_(errs);
}

QtCharts::QChartView* MainWindow::createChartInPlaceholder(QWidget* ph, QtCharts::QChart* chart)
{
  using namespace QtCharts;
  auto* view = new QChartView(chart, ph);
  view->setRenderHint(QPainter::Antialiasing);
  view->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  auto* lay = new QVBoxLayout(ph);
  lay->setContentsMargins(0,0,0,0);
  lay->addWidget(view);
  return view;
}

void MainWindow::setupCalibrationCharts_()
{
  using namespace QtCharts;

  // ---- Smile IV(K) ----
  {
    smileChart_ = new QChart();
    smileChart_->setTitle("Smile IV(K)");
    smileChart_->setMargins(QMargins(8,8,18,12));
    smileChart_->legend()->setVisible(true);

    // axes
    axSmileX_ = new QValueAxis(); axSmileX_->setTitleText("Strike K"); axSmileX_->setLabelFormat("%.0f");
    axSmileY_ = new QValueAxis(); axSmileY_->setTitleText("IV");       axSmileY_->setLabelFormat("%.3f");
    smileChart_->addAxis(axSmileX_, Qt::AlignBottom);
    smileChart_->addAxis(axSmileY_, Qt::AlignLeft);

    // séries
    smileLine_  = new QLineSeries();     smileLine_->setName("spline");
    smileNodes_ = new QScatterSeries();  smileNodes_->setName("nodes"); smileNodes_->setMarkerSize(6.0);

    smileChart_->addSeries(smileLine_);
    smileChart_->addSeries(smileNodes_);
    smileLine_->attachAxis(axSmileX_);  smileLine_->attachAxis(axSmileY_);
    smileNodes_->attachAxis(axSmileX_); smileNodes_->attachAxis(axSmileY_);

    chartCalibSmile_ = createChartInPlaceholder(ui->chartCalibSmilePH, smileChart_);
  }

  // ---- Résidus ΔP ----
  {
    residChart_ = new QChart();
    residChart_->setTitle("Résidus (price_mkt - price_fit)");
    residChart_->setMargins(QMargins(8,8,18,12));
    residChart_->legend()->setVisible(true);

    // axes persistants
    axResidX_ = new QBarCategoryAxis(); axResidX_->setTitleText("Quote #");
    axResidY_ = new QValueAxis();       axResidY_->setTitleText("Erreur prix"); axResidY_->setLabelFormat("%.3f");
    residChart_->addAxis(axResidX_, Qt::AlignBottom);
    residChart_->addAxis(axResidY_, Qt::AlignLeft);

    // série (vide au départ)
    residBars_ = new QBarSeries();
    residChart_->addSeries(residBars_);
    residBars_->attachAxis(axResidX_);
    residBars_->attachAxis(axResidY_);

    chartCalibResid_ = createChartInPlaceholder(ui->chartCalibResidPH, residChart_);
  }
}

void MainWindow::updateSmileChartForT_(double T)
{
  if (!smileChart_) return;

  QVector<QPointF> pts;
  double kmin=1e300, kmax=-1e300, ivmin=1e300, ivmax=-1e300;

  for (const auto& r : msCalibCache_.rows) {
    if (std::fabs(r.T - T) < 1e-10 && std::isfinite(r.iv_mid) && r.iv_mid > 0.0) {
      pts.push_back(QPointF(r.K, r.iv_mid));              // <-- QPointF(...)
      kmin = std::min(kmin, r.K); kmax = std::max(kmax, r.K);
      ivmin= std::min(ivmin, r.iv_mid); ivmax= std::max(ivmax, r.iv_mid);
    }
  }

  QVector<QPointF> curve;
  if (kmin < kmax) {
    const int NK = 200;
    for (int i=0; i<NK; ++i) {
      const double u = double(i)/(NK-1);
      const double K = std::exp(std::log(kmin) + u*(std::log(kmax)-std::log(kmin)));
      const double iv = smileCalibCache_.iv_at(T, K);
      if (std::isfinite(iv)) {
        curve.push_back(QPointF(K, iv));                  // <-- QPointF(...)
        ivmin = std::min(ivmin, iv); ivmax = std::max(ivmax, iv);
      }
    }
  }

  // noms cohérents avec tes membres
  smileNodes_->replace(pts);                               // <-- smileNodes_ (pas smilePts_)
  smileLine_->replace(curve);

  if (kmin < kmax) {
    if (axSmileX_) axSmileX_->setRange(kmin, kmax);
    if (axSmileY_) axSmileY_->setRange(std::max(0.0, ivmin*0.9), ivmax*1.1);

  }
}

// --- Crée/initialise les charts du tab Calibration (placeholders du .ui)
void MainWindow::setupCalibrationTab_()
{
  // Tu utilises des placeholders du .ui -> on ne crée pas l’onglet,
  // on instancie seulement les QChart dans ces containers.
  setupCalibrationCharts_();
}

// --- Slot UI : bouton "Construire Smile" -> worker
void MainWindow::onCalibBuildSmile()
{
  startCalibStage_("Construction du smile…", 0);
  QMetaObject::invokeMethod(calib_, "buildSmile", Qt::QueuedConnection);
}

void MainWindow::updateResidChart_(const std::vector<double>& errs)
{
  residChart_->removeSeries(residBars_);
  residBars_ = new QtCharts::QBarSeries(residChart_);

  auto* set = new QtCharts::QBarSet("ΔP");
  QStringList cats;
  double emin=0.0, emax=0.0;
  for (size_t i=0;i<errs.size();++i) {
    *set << errs[i];
    cats << QString::number(i+1);
    emin = std::min(emin, errs[i]);
    emax = std::max(emax, errs[i]);
  }

  residBars_->append(set);
  residChart_->addSeries(residBars_);
  residBars_->attachAxis(axResidX_);
  residBars_->attachAxis(axResidY_);

  if (axResidX_) {
  QStringList cats; cats.reserve((int)errs.size());
  for (int i=0;i<(int)errs.size();++i) cats << QString::number(i+1);
  axResidX_->clear(); axResidX_->append(cats);
  }
  if (axResidY_) {
    double mn=0, mx=0;
    for (double e: errs){ mn = std::min(mn,e); mx = std::max(mx,e); }
    const double pad = 0.05 * std::max(1.0, std::fabs(mx - mn));
    axResidY_->setRange(mn - pad, mx + pad);
  }
}

static QString safePercentText_(int cur, int tot) {
  if (tot <= 0) return QString();
  const int pct = int(std::round(100.0 * double(std::max(0,cur)) / double(std::max(1,tot))));
  return QString("%1%").arg(pct);
}

QString MainWindow::fmtEta_(double s) {
  if (s < 1.0) return QString::number(std::max(0.1, s), 'f', 1) + "s";
  int sec = int(std::round(s));
  int m = sec / 60; sec %= 60;
  if (m == 0) return QString("%1s").arg(sec);
  return QString("%1m%2").arg(m).arg(sec, 2, 10, QLatin1Char('0'));
}

void MainWindow::setCalibProgressStyle_(CalibProgState st) {
  calibProgState_ = st;
  auto* p = ui->calibProgress;
  switch (st) {
    case CalibProgState::Busy:
      p->setStyleSheet("QProgressBar{border:1px solid #bbb;border-radius:3px;text-align:center;}"
                       "QProgressBar::chunk{background:#E53935;}"); // rouge
      break;
    case CalibProgState::Done:
      p->setStyleSheet("QProgressBar{border:1px solid #bbb;border-radius:3px;text-align:center;}"
                       "QProgressBar::chunk{background:#43A047;}"); // vert
      break;
    case CalibProgState::Error:
      p->setStyleSheet("QProgressBar{border:1px solid #bbb;border-radius:3px;text-align:center;}"
                       "QProgressBar::chunk{background:#FB8C00;}"); // orange
      break;
    case CalibProgState::Idle:
    default:
      p->setStyleSheet("");
      break;
  }
}

void MainWindow::startCalibStage_(const QString& label, int total) {
  ui->lblCalibMeta->setText(label);
  calibEtaTimer_.restart();
  ui->calibProgress->setFormat(label + "  %p%");
  if (total > 0) { ui->calibProgress->setRange(0, total); ui->calibProgress->setValue(0); }
  else           { ui->calibProgress->setRange(0, 0); } // indéterminé
  setCalibProgressStyle_(CalibProgState::Busy);
}

void MainWindow::updateCalibStage_(int cur, int total) {
  auto* p = ui->calibProgress;
  if (total > 0) {
    p->setRange(0, total);
    p->setValue(cur);
    // ETA = elapsed * (total/cur - 1)
    const qint64 ms = calibEtaTimer_.elapsed();
    if (cur > 0) {
      const double elapsed = ms / 1000.0;
      const double eta = elapsed * (double(total) / double(cur) - 1.0);
      const QString txt = QString("%1  %2  ETA %3")
            .arg(ui->lblCalibMeta->text(),
                 safePercentText_(cur,total),
                 fmtEta_(std::max(0.0, eta)));
      p->setFormat(txt);
    }
  } else {
    // indéterminé: juste garder le label en rouge
    p->setRange(0,0);
    p->setFormat(ui->lblCalibMeta->text());
  }
}

void MainWindow::finishCalibStage_(bool ok, const QString& label) {
  if (!label.isEmpty()) ui->lblCalibMeta->setText(label);
  ui->calibProgress->setRange(0, 1);
  ui->calibProgress->setValue(1);
  ui->calibProgress->setFormat(label.isEmpty() ? "%p%" : (label + "  %p%"));
  setCalibProgressStyle_(ok ? CalibProgState::Done : CalibProgState::Error);
}

void MainWindow::showFileToast_(const QString& label, const QString& path)
{
  // Nettoie un éventuel widget précédent
  // (on peut aussi garder plusieurs toasts si tu préfères)
  statusBar()->clearMessage();

  auto* wrap = new QWidget(this);
  auto* lay  = new QHBoxLayout(wrap);
  lay->setContentsMargins(6,0,6,0);

  auto* info = new QLabel(QString("%1: %2")
                          .arg(label,
                               QDir::toNativeSeparators(path)), wrap);
  auto* btnOpen   = new QPushButton("Ouvrir", wrap);
  auto* btnFolder = new QPushButton("Ouvrir dossier", wrap);
  auto* btnCopy   = new QPushButton("Copier", wrap);

  lay->addWidget(info);
  lay->addSpacing(8);
  lay->addWidget(btnOpen);
  lay->addWidget(btnFolder);
  lay->addWidget(btnCopy);

  statusBar()->addPermanentWidget(wrap, 1);

  // Actions
  connect(btnOpen, &QPushButton::clicked, this, [path]{
    QDesktopServices::openUrl(QUrl::fromLocalFile(path));
  });
  connect(btnFolder, &QPushButton::clicked, this, [path]{
    const QString dir = QFileInfo(path).absolutePath();
    QDesktopServices::openUrl(QUrl::fromLocalFile(dir));
  });
  connect(btnCopy, &QPushButton::clicked, this, [this, path]{
    QGuiApplication::clipboard()->setText(path);
    statusBar()->showMessage("Chemin copié dans le presse-papiers", 2000);
  });

  // Auto-disparition après 10s (facultatif)
  QTimer::singleShot(10000, this, [this, wrap]{
    statusBar()->removeWidget(wrap);
    wrap->deleteLater();
  });
}

QJsonObject MainWindow::nodeToJson(double K, double iv) {
  return QJsonObject{{"K", K}, {"iv", iv}};
}

QJsonArray MainWindow::serializeSmileSlices_() const {
  QJsonArray out;
  if (smileCalibCache_.slices.empty()) return out;

  // Liste des maturités de la surface
  QVector<double> Ts;
  Ts.reserve((int)smileCalibCache_.slices.size());
  for (const auto& sl : smileCalibCache_.slices) Ts.push_back(sl.T);

  for (double T : Ts) {
    double kmin =  std::numeric_limits<double>::infinity();
    double kmax = -std::numeric_limits<double>::infinity();

    // borne Kmin/Kmax à cette T à partir du marché
    for (const auto& r : msCalibCache_.rows) {
      if (std::fabs(r.T - T) < 1e-10 && std::isfinite(r.K) && r.K > 0.0) {
        kmin = std::min(kmin, r.K);
        kmax = std::max(kmax, r.K);
      }
    }
    // fallback: toute la table si aucune quote pile à T
    if (!(kmin < kmax)) {
      for (const auto& r : msCalibCache_.rows) {
        if (std::isfinite(r.K) && r.K > 0.0) {
          kmin = std::min(kmin, r.K);
          kmax = std::max(kmax, r.K);
        }
      }
    }
    if (!(kmin < kmax)) continue;

    // Échantillonnage de la courbe IV(K)
    constexpr int NK = 25;
    QJsonArray nodes;
    const double logKmin = std::log(kmin);
    const double logKmax = std::log(kmax);

    for (int i = 0; i < NK; ++i) {
      const double u  = (NK == 1 ? 0.0 : double(i) / double(NK - 1));
      const double K  = std::exp(logKmin + u * (logKmax - logKmin));
      const double iv = smileCalibCache_.iv_at(T, K);
      if (std::isfinite(iv) && iv > 0.0) {
        nodes.append(QJsonObject{{"K", K}, {"iv", iv}});
      }
    }

    out.append(QJsonObject{
      {"T", T},
      {"nodes", nodes}
    });
  }
  return out;
}

bool MainWindow::saveWidgetPng(QWidget* w, const QString& outPath,
                               const QSize& targetPx, qreal dpr) {
  if (!w) return false;
  const QSize base = targetPx.isValid() ? targetPx : w->size();
  const QSize hi(int(base.width() * dpr), int(base.height() * dpr));

  QImage img(hi, QImage::Format_ARGB32_Premultiplied);
  img.setDevicePixelRatio(dpr);
  img.fill(Qt::transparent);

  w->render(&img);           // rend le widget dans l'image
  return img.save(outPath, "PNG");
}

void MainWindow::onExportChartsPng() {
  // 1) choisir un dossier parent
  const QString parent = QFileDialog::getExistingDirectory(
      this, tr("Choisir un dossier parent d'export"));
  if (parent.isEmpty()) return;

  // 2) construire un sous-dossier unique (horodaté + nom projet si dispo)
  auto sanitize = [](QString s){
    s.replace(QRegExp(R"([^A-Za-z0-9_\-\.])"), "_");
    return s.left(40);
  };
  const QString ts  = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
  const QString proj = projectLabel_ ? sanitize(projectLabel_->text()) : "RW";
  const QString sub  = QString("%1_%2").arg(ts, proj.isEmpty() ? "RW" : proj);
  QDir parentDir(parent);
  const QString outDirPath = parentDir.filePath(sub);

  QDir outDir;
  if (!outDir.mkpath(outDirPath)) {
    statusBar()->showMessage(tr("Impossible de créer le dossier: %1").arg(outDirPath), 5000);
    return;
  }

  // 3) lister les widgets à exporter (adapte les noms si besoin)
  struct Item { QWidget* w; const char* base; };
  const Item items[] = {
    { ui->chartContainer,   "pricing_convergence" },
    { ui->greeksChartContainer, "greeks" },
    { ui->stressChartContainer, "stress" },
    { ui->chartCalibSmilePH,  "calib_smile" },
    { ui->chartCalibResidPH,  "calib_residuals" },
    { ui->lblCalibHeatmap,    "calib_heatmap" }
  };

  // 4) exporter
  int nOK = 0, nAll = 0;
  for (const Item& it : items) {
    if (!it.w) continue;
    ++nAll;
    const QString pngPath = QDir(outDirPath).filePath(QString("%1.png").arg(it.base));
    if (saveWidgetPng(it.w, pngPath, /*targetPx=*/QSize(), /*dpr=*/2.0)) ++nOK;
  }

  statusBar()->showMessage(
      tr("Export PNG: %1/%2 fichiers dans %3").arg(nOK).arg(nAll).arg(outDirPath), 6000);
}
