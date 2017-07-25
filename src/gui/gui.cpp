#include "mainwindow.h"
#include <QApplication>

#ifdef ESPRESOGUIPROJECTBUILD
#include "../config/ecf/environment.h"
#endif

int main(int argc, char *argv[])
{
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Init(&argc, &argv);
#endif
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    a.exec();
#ifdef ESPRESOGUIPROJECTBUILD
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
}
