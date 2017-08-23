#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QStandardItemModel>

#include "../config/ecf/ecf.h"
#include <iostream>

using namespace espreso;

template <typename Ttype>
std::ostream& operator<<(std::ostream& os, const std::vector<Ttype> &v)
{
	for (size_t i = 0; i < v.size(); i++) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

static void printECF(const ECFObject &object, size_t indent)
{
	auto printindent = [&] () {
		for (size_t j = 0; j < indent; j++) {
			std::cout << " ";
		}
	};

	auto printparameter = [&] (const ECFParameter *parameter) {
		if (parameter == NULL || !parameter->metadata.isallowed()) {
			return;
		}
		printindent();
		std::cout << parameter->name << " ";
		if (parameter->isValue()) {
			std::cout << " = " << parameter->getValue() << ";\n";
			if (parameter->metadata.datatype.front() == ECFDataType::EXPRESSION) {
				printindent();
				std::cout << "EXPRESSION: " << parameter->metadata.variables;
			}
			if (parameter->metadata.tensor != NULL) {
				printindent();
				std::cout << "TENSOR: " << parameter->metadata.tensor->size << "x" << parameter->metadata.tensor->size << "\n";
			}
		} else if (parameter->isObject()) {
			std::cout << "{\n";
			printECF(*dynamic_cast<const ECFObject*>(parameter), indent + 2);
			printindent();
			std::cout << "}\n";
		} else {
			// Separators etc..
			std::cout << "\n";
		}

	};

	for (size_t i = 0; i < object.parameters.size(); i++) {
		printparameter(object.parameters[i]);
	}
};

MainWindow::MainWindow(int *argc, char ***argv)
: QMainWindow(0), ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	espreso::ECFConfiguration ecf(argc, argv);
	printECF(ecf, 0);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_pushButton_clicked()
{
	ui->treeView->move(QPoint(20, 20));
}

void MainWindow::mousePressEvent(QMouseEvent *ev)
{
	last = ev->pos();
}

void MainWindow::mouseMoveEvent(QMouseEvent *ev)
{
	QPoint diff = last - ev->pos();
	ui->treeView->move(ui->treeView->pos() - diff);
	last = ev->pos();
}
