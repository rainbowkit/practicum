#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <iostream>
#include <cmath>

#include "task.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow

{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
private slots:
    void calculate();
    void replot();
private:
    void rescaleAndMatch();

    Ui::MainWindow *ui;
    QVector<double> x, y;
};

#endif // MAINWINDOW_H
