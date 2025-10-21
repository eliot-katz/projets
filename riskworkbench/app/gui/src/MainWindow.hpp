#pragma once
#include <QMainWindow>
#include <QThread>
#include <memory>
#include <optional>
#include <cstddef>
#include <vector>

#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QBarCategoryAxis>
#include <QtCharts/QCategoryAxis>
#include <QtCharts/QStackedBarSeries>

#include <QTimer>
#include <QLabel>
#include <QShortcut>
#include <QJsonObject>

#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/calibration/calib.hpp>
#include <rw/market/surface.hpp>  
#include <rw/smile/smile.hpp>     
#include <rw/qc/qc.hpp> 


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

namespace gui { class McWorker; }
namespace gui { class CalibWorker; }

namespace QtCharts {
  class QChartView;
  class QChart;
  class QLineSeries;
  class QAreaSeries;
  class QValueAxis;
}

class QTableWidget;
class QProgressBar;
class QComboBox;
class QLabel;
class QPushButton;

class MainWindow : public QMainWindow {
  Q_OBJECT
public:
  explicit MainWindow(QWidget* parent = nullptr);
  ~MainWindow() override;   // destructeur

private slots:
  // Tes handlers "plats"
  void onRunMC();
  void onRunGreeks();
  void onRunStress();
  void onLoadDefaults();
  void onReset();

  // Callbacks worker MC
  void onMcProgress(std::size_t n, double mean, double half);
  void onMcFinished(double price, double se, double lo, double hi,
                    std::size_t nEff, long long ms);
  void onMcFailed(const QString& why);
  void onMcCanceled();

  void onSnapBaseline(); 

  void onSaveProject();
  void onLoadProject();

  void onShortcutSave();    // Ctrl+S
  void onShortcutSaveAs();  // Ctrl+Shift+S

  // --- Calibration UI ---
  // actions UI
  void onCalibOpenCsv();
  void onCalibInvertIv();
  void onCalibBuildSmile();
  void onCalibFitGlobal();
  void onCalibRunQc();
  void onCalibExportCsv();
  void onCalibSaveProject();
  void onCalibSelectT(int idx);

  // retours worker
  void onCalibMessage(const QString& m);
  void onCalibProgress(const QString& stage, int cur, int tot);
  void onCalibLoaded(size_t nRows, double S0, double r, double q, const QString& und);
  void onCalibInverted(size_t priced, size_t ok, size_t bad);
  void onCalibSmileBuilt(const QVector<double>& maturities);
  void onCalibGlobalFitted(const rw::calib::CalibReport& rep, double sigma);
  void onCalibQcDone(double rmse_price, double rmse_iv, size_t n,
                   const std::vector<rw::qc::ErrorRow>& rows);

  void onCalibExported(const QString& path);
  void onCalibProjectSaved(const QString& path);
  void onCalibFailed(const QString& why);

  void onExportChartsPng();


private:
  Ui::MainWindow* ui;

  // Worker thread pour MC
  QThread* mcThread_ {nullptr};
  gui::McWorker* mcWorker_ {nullptr};

  // Wiring des signaux UI
  void wireSignals();

  // Helpers UI existants
  void setPricingResults(double price, double se, double ciLow, double ciHigh,
                         std::size_t nEff, long long ms);
  void setGreeksResults(double d, double v, double g, double r, double t);

  // Gestion worker
  void startMcWorker();
  void stopMcWorker(bool hard = false);

  // Utilitaires
  void clearConvergenceTable();
  void appendConvergenceRow(std::size_t n, double est, double half);

  // ===== Convergence chart =====
  void setupConvergenceChart();       // create chart in the placeholder
  void resetConvergenceChart();       // clear series at start of a run
  void setupConvergenceLegend();
  void updateAxesForPoint(qreal x, qreal yLow, qreal yHigh);

  QtCharts::QChartView*  convChartView_{nullptr};
  QtCharts::QLineSeries* meanSeries_{nullptr};
  QtCharts::QLineSeries* ciUpper_{nullptr};
  QtCharts::QLineSeries* ciLower_{nullptr};
  QtCharts::QAreaSeries* ciBand_{nullptr};
  QtCharts::QValueAxis*  xAxis_{nullptr};
  QtCharts::QValueAxis*  yAxis_{nullptr};

  // Lignes/séries pour le graphe de convergence
  QtCharts::QLineSeries*    bsLine_  {nullptr};  // ligne de référence Black–Scholes
  QtCharts::QScatterSeries* lastPt_  {nullptr};  // marqueur du dernier point

  double yMinSeen_ =  std::numeric_limits<double>::infinity();
  double yMaxSeen_ = -std::numeric_limits<double>::infinity();


  // Baseline pour le stress
  std::optional<rw::market::MarketData> baseMkt_;
  std::optional<rw::market::Instrument> baseInst_;
  double                  baseSigma_{0.20};
  bool                    baselineSet_{false};

  QJsonObject makeProjectJson() const;
  void        loadProjectJson(const QJsonObject& obj);
  QString projectsDir() const;  // retourne le chemin absolu vers data/projects

  // --- Projet courant ---
  QLabel*  projectLabel_{nullptr};     // dans la status bar
  QString  currentProjectPath_;        // chemin complet du projet
  bool     projectDirty_{false};       // indicateur "modifié"

  void setCurrentProject(const QString& path, bool dirty=false);
  void markProjectDirty();             // à appeler quand l'UI change
  QString projectDisplayName() const;  // "project.json" ou "No project"

  // --- Raccourcis & helpers Save ---
  void setupProjectShortcuts();
  bool writeProjectTo(const QString& path, QString* errMsg=nullptr) const;

  // ===== Greeks bar chart =====
  QtCharts::QChartView*      greeksChartView_{nullptr};
  QtCharts::QChart*          greeksChart_{nullptr};
  QtCharts::QBarSeries*      greeksSeries_{nullptr};
  QtCharts::QBarSet*         greeksSet_{nullptr};        // 5 valeurs: Δ, Vega, Γ, ρ, Θ
  QtCharts::QBarCategoryAxis* gAxisX_{nullptr};
  QtCharts::QValueAxis*       gAxisY_{nullptr};

  void setupGreeksChart();
  void resetGreeksChart();
  void updateGreeksChart(double d,double v,double g,double r,double t);
  double greeksMinSeen_{0.0};
  double greeksMaxSeen_{0.0};

  // Auto-run (debounce)
  QTimer* stressDebounce_{nullptr};
  bool    stressBusy_{false};
  bool    stressPending_{false};
  void    armStressDebounce_(int ms = 300);

  // Chart Stress: barres “MC vs Taylor” + breakdown
  QtCharts::QChartView*   stressChartView_{nullptr};
  QtCharts::QChart*       stressChart_{nullptr};
  QtCharts::QStackedBarSeries*   stressTaylorSeries_{nullptr}; // stacked (Δ, Γ, Vega, Rho)
  QtCharts::QBarSeries*   stressMCSeries_{nullptr};     // barre unique MC
  QtCharts::QValueAxis*   stressAxisY_{nullptr};
  QtCharts::QBarCategoryAxis* stressAxisX_{nullptr};

  // Pour construire/rafraîchir
  void initStressControls();     // sliders↔spins + auto-run
  void setupStressChart();       // crée le chart
  void resetStressChart();       // vide séries
  void setupStressLegend();      // construit la légende custom
  void updateStressChart(double pnlMC,
                         double d_contrib, double g_contrib,
                         double v_contrib, double r_contrib);

  // --- Validation/Test harness ---
  enum class TestKind { PricingBS, Convergence, GreeksBS, StressCRN, SaveReload };

  struct TestRow {
    QString test, metric, value, threshold, result; // pour la table
    bool pass{false};
  };

  void runAllTests();
  void runTestPricingBS();      // BS ∈ [CI]
  void runTestConvergence();    // pente log(IC) ~ -0.5
  void runTestGreeksBS();       // z-score ≤ 2 (voir recette)
  void runTestStressCRN();      // Var(P&L) << Var(price) (voir recette)
  void runTestSaveReload();     // Idempotence JSON

  void testsClearTable() const;
  void testsAppendRow(const TestRow& r) const;
  void testsLog(const QString& line) const;

  // Convergence accumulation
  std::vector<double> convN_;
  std::vector<double> convHalf_;

  // Petit utilitaire régression y = a + b x (retourne b)
  static double slopeLeastSquares(const std::vector<double>& x, const std::vector<double>& y);
  static double sampleVariance(const std::vector<double>& x, bool unbiased = true);

  /////////////  Calibration   /////////////
  void setupCalibrationTab_();
  void wireCalibration_();
  void repaintSmile_();
  void repaintResiduals_();
  void repaintHeatmap_();

  // ===== Calibration charts =====
  QtCharts::QChartView*  chartCalibSmile_{nullptr};
  QtCharts::QChart*      smileChart_{nullptr};
  QtCharts::QLineSeries* smileLine_{nullptr};
  QtCharts::QScatterSeries* smileNodes_{nullptr};

  QtCharts::QChartView*  chartCalibResid_{nullptr};
  QtCharts::QChart*      residChart_{nullptr};
  QtCharts::QBarSeries*  residBars_{nullptr};

  // Smile
  QtCharts::QValueAxis* axSmileX_{nullptr};
  QtCharts::QValueAxis* axSmileY_{nullptr};

  // Résidus
  QtCharts::QBarCategoryAxis* axResidX_{nullptr};
  QtCharts::QValueAxis*       axResidY_{nullptr};


  // helpers
  QtCharts::QChartView* createChartInPlaceholder(QWidget* ph, QtCharts::QChart* chart);
  void setupCalibrationCharts_();
  void updateSmileChartForT_(double T);                         // lit msCalibCache_/smileCalibCache_
  void updateResidChart_(const std::vector<double>& residuals); // valeurs ΔP (price_mkt - price_fit)


  // worker
  QThread* calibThread_{nullptr};
  gui::CalibWorker* calib_{nullptr};

  // caches
  rw::market::MarketSurface msCalibCache_;
  rw::smile::SmileSurface   smileCalibCache_;
  std::vector<rw::qc::ErrorRow> residCalibCache_;
  QVector<double> maturitiesCalib_;

  // --- Calibration: progress/ETA ---
  enum class CalibProgState { Idle, Busy, Done, Error };

  QElapsedTimer calibEtaTimer_;
  CalibProgState calibProgState_{CalibProgState::Idle};

  void setCalibProgressStyle_(CalibProgState st);
  void startCalibStage_(const QString& label, int total);     // reset timer + rouge
  void updateCalibStage_(int cur, int total);                  // met % + ETA
  void finishCalibStage_(bool ok=true, const QString& label="");
  static QString fmtEta_(double seconds);                      // "3.2s", "1m05"

  void showFileToast_(const QString& label, const QString& path);
  void applySmileSettingsToWorker_();

  static QJsonObject nodeToJson(double K, double iv);
  QJsonArray serializeSmileSlices_() const; 

  static bool saveWidgetPng(QWidget* w,
                            const QString& outPath,
                            const QSize& targetPx = QSize(),
                            qreal devicePixelRatio = 2.0);
  
};
