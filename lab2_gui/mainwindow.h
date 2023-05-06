#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <iostream>
#include <cmath>

#include "qmessagebox.h"
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
    void inverseOnToggle();
private:
    void rescaleAndMatch();
    void alert(const QString prompt, QMessageBox::Icon icon);

    bool inverse;
    QMessageBox msgBox;
    Ui::MainWindow *ui;
    QVector<double> x, y;
};

#endif // MAINWINDOW_H
