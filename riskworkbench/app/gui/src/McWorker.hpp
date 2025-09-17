#pragma once
#include <QObject>
#include <atomic>
#include <cstddef>
#include <QString>

// RW types (on passe par valeur ⇒ on inclut ici)
#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/models/gbm.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/config/greeks_config.hpp>

namespace gui {

class McWorker : public QObject {
  Q_OBJECT
public:
  explicit McWorker(QObject* parent = nullptr);
  ~McWorker() override = default;

public slots:
  // Exécution Monte Carlo (mono-thread), avec progress() par batch.
  void runPricing(rw::market::MarketData mkt,
                  rw::market::Instrument inst,
                  rw::config::McConfig cfg,
                  double sigma);

  // Lance un calcul de greeks selon "method" = "brv" | "pw" | "lrm"
  void runGreeks(rw::market::MarketData mkt,
                 rw::market::Instrument inst,
                 rw::config::McConfig mc,
                 rw::config::GreeksConfig gk,
                 double sigma,
                 QString method);

  // run stress — baseline & one stressed scenario (both repriced by MC)
  void runStress(rw::market::MarketData baseMkt,
                 rw::market::Instrument baseInst,
                 rw::config::McConfig cfg,
                 double baseSigma,
                 double dS_rel,   // relative shock on S0 (e.g. +0.05 => +5%)
                 double dSigma,   // absolute shock on sigma (e.g. +0.02)
                 double dR,       // absolute shock on r (e.g. +0.001)
                 double dT);      // absolute shock on T in years (can be +/-)

  // Demande d’arrêt asynchrone
  void requestStop();

signals:
  void progress(std::size_t nDone, double mean, double halfwidth95);
  void finished(double price, double std_error,
                double ci_low, double ci_high,
                std::size_t n_effective, long long elapsed_ms);
  void failed(QString why);
  void canceled();
  void greeksFinished(double delta, double vega, double gamma,
                      double rho, double theta, long long elapsed_ms);
  // stress results
  void stressFinished(double basePrice, double stressedPrice,
                      double pnl, long long elapsed_ms);

private:
  std::atomic<bool> stop_{false};
};

} // namespace gui
