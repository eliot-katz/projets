#pragma once
#include <QtCharts/QChartView>
namespace QtCharts { class QChartView; }
class ChartView : public QtCharts::QChartView {
public:
  explicit ChartView(QWidget* parent=nullptr) : QtCharts::QChartView(parent) {}
};
