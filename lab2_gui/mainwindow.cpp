
#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Setup data for plot
    x = QVector<double>(101);
    y = QVector<double>(101);

    // Information on plot
    ui->graph->legend->setVisible(true);

    // Main graph and dots
    ui->graph->addGraph();
    ui->graph->addGraph();
    ui->graph->addGraph();
    ui->graph->addGraph();

    ui->graph->graph(0)->setPen(QColor(255, 0, 0, 100));
    ui->graph->graph(0)->setName("f(x) = sin(x) – x^2/2");

    // Make graph look like dots
    ui->graph->graph(1)->setPen(QColor(0, 150, 255, 255));
    ui->graph->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->graph->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    ui->graph->graph(1)->setName("Ближайшие узлы");

    ui->graph->graph(2)->setPen(QColor(0, 0, 0, 255));
    ui->graph->graph(2)->setLineStyle(QCPGraph::lsNone);
    ui->graph->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->graph->graph(2)->setName("Остальные узлы");

    // Point of interpolation
    ui->graph->graph(3)->setPen(QColor(255, 0, 0, 255));
    ui->graph->graph(3)->setLineStyle(QCPGraph::lsNone);
    ui->graph->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 6));
    ui->graph->graph(3)->setName("Точка интерпол.");

    // Additional properties
    ui->graph->xAxis->setLabel("x");
    ui->graph->yAxis->setLabel("y");
    ui->graph->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    // Layout
    ui->graph->xAxis->setRange(ui->rightBound->value(), ui->leftBound->value());
    replot();
    rescaleAndMatch();
}

void MainWindow::rescaleAndMatch() {
    ui->graph->graph(0)->rescaleValueAxis();
    double currentStart = ui->graph->xAxis->range().lower;
    ui->graph->yAxis->setRange(currentStart, currentStart + ui->graph->xAxis->range().size());
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::calculate() {
    const double x0 = ui->xInter->value();

    // Values validation
    const int nodesCount = ui->nodesCount->value();
    const int power = ui->power->value();
    if (power >= nodesCount) {
        ui->power->setStyleSheet("border: 1px solid red;");
        return;
    } else {
        ui->power->setStyleSheet("");
    }

    const double a = ui->leftBound->value();
    const double b = ui->rightBound->value();
    if (a >= b) {
        ui->rightBound->setStyleSheet("border: 1px solid red;");
        ui->leftBound->setStyleSheet("border: 1px solid red;");
        return;
    } else {
        ui->rightBound->setStyleSheet("");
        ui->leftBound->setStyleSheet("");
    }

    // Dots
    auto table = generateTable(a, b, nodesCount);
    sort(table.begin(), table.end(), CloserTo(x0));
    QVector<double> x1, y1;
    ui->nodes->clear();
    for (int i = 0; i < power + 1; ++i) {
        x1.push_back(table.at(i).first);
        y1.push_back(table.at(i).second);
        ui->nodes->addItem(QString::fromStdString("X =  ") +
                           QString::number(table.at(i).first) +
                           QString::fromStdString("; Y =  ") +
                           QString::number(table.at(i).second));
    }
    ui->graph->graph(1)->setData(x1, y1);
    x1.clear();
    y1.clear();
    ui->nodes->addItem(QString::fromStdString(" "));
    for (int i = power + 1; i < table.size(); ++i) {
        x1.push_back(table.at(i).first);
        y1.push_back(table.at(i).second);
        ui->nodes->addItem(QString::fromStdString("X =  ") +
                           QString::number(table.at(i).first) +
                           QString::fromStdString("; Y =  ") +
                           QString::number(table.at(i).second));
    }
    ui->graph->graph(2)->setData(x1, y1);
    x1.clear();
    y1.clear();

    // Interpolation itself
    auto result1 = lagrange(x0, table, power);
    auto result2 = newton(x0, table, power);
    ui->result->setText(QString::fromStdString("Форма Лагранжа: ") + QString::number(result1.first) +
                        QString::fromStdString("\nПогрешность: ") + QString::number(result1.second) +
                        QString::fromStdString("\nФорма Ньютона: ") + QString::number(result2.first) +
                        QString::fromStdString("\nПогрешность: ") + QString::number(result2.second));

    // Plotting result in Lagrange form
    x1.push_back(x0);
    y1.push_back(result1.first);
    ui->graph->graph(3)->setData(x1, y1);

    replot();
}

void MainWindow::replot() {
    auto axis = ui->graph->xAxis->range();

    // Main graph
    for (int i=0; i<101; ++i)
    {
        x[i] = axis.lower + i/100.0 * axis.size();
        y[i] = sin(x[i]) - x[i]*x[i]/2;
    }
    ui->graph->graph(0)->setData(x, y);

    ui->graph->replot();
}
