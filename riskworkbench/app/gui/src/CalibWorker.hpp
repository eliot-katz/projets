#pragma once
#include <QObject>
#include <QString>
#include <QVector>
#include <vector>

#include <rw/market/surface.hpp>
#include <rw/calibration/calib.hpp>
#include <rw/smile/smile.hpp>
#include <rw/qc/qc.hpp>


namespace gui {

class CalibWorker : public QObject {
  Q_OBJECT
public:
  explicit CalibWorker(QObject* parent=nullptr);

public slots:
  void loadCsv(const QString& path);
  void invertIv();                       // complète iv_mid via inversion BS
  void buildSmile();                     // lissage slices → SmileSurface
  void fitGlobal(double sigma0=0.2);     // σ global
  void runQc(double tick=0.0, double onlyConsistentTol=-1.0, bool priceFromIv=false);
  void exportSurfaceCsv(const QString& outDirBase = QString("data/iv_exports"));
  void saveProjectJson();

signals:
  void message(const QString& text);
  void progress(const QString& stage, int cur, int total);

  void loadedMarket(quint64 nRows, double S0, double r, double q, const QString& underlying);
  void inverted(quint64 priced, quint64 ok, quint64 bad);
  void qcDone(double rmse_price, double rmse_iv, quint64 n, const std::vector<rw::qc::ErrorRow>& rows);
  void smileBuilt(const QVector<double>& maturities);
  void globalFitted(const rw::calib::CalibReport& rep, double sigma);
  

  void exportedCsv(const QString& path);
  void projectSaved(const QString& path);

  void failed(const QString& why);

public:
  const rw::market::MarketSurface& market() const { return ms_; }
  const rw::smile::SmileSurface&   smile()  const { return smile_; }
  const std::vector<rw::qc::ErrorRow>& residuals() const { return residuals_; }
  QString lastUnderlying() const { return underlying_; }

private:
  rw::market::MarketSurface ms_;
  rw::smile::SmileSurface   smile_;
  std::vector<rw::qc::ErrorRow> residuals_;
  QString underlying_{"UNKNOWN"};

  static QString sanitize(const QString&);

};

} // namespace gui
