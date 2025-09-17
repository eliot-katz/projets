#pragma once
#include <QMainWindow>
#include <QThread>
#include <memory>
#include <optional>
#include <cstddef>

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


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

namespace gui { class McWorker; }

namespace QtCharts {
  class QChartView;
  class QChart;
  class QLineSeries;
  class QAreaSeries;
  class QValueAxis;
}

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


  
};
