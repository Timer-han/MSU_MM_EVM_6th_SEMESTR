#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong input arguments!", "Wrong input arguments!");
        return -1;
    }

    action = tool_bar->addAction("Change &Function", graph_area, SLOT(change_f()));
    action->setShortcut(QString("0"));
    
    action = tool_bar->addAction("Change Show &Mode", graph_area, SLOT(change_show_mode()));
    action->setShortcut(QString("1"));
    
    action = tool_bar->addAction("&Compress Area", graph_area, SLOT(decrease_visible_area()));
    action->setShortcut(QString("2"));
    
    action = tool_bar->addAction("&Expand Area", graph_area, SLOT(increase_visible_area()));
    action->setShortcut(QString("3"));
    
    action = tool_bar->addAction("n+", graph_area, SLOT(increase_triangulation()));
    action->setShortcut(QString("4"));
    
    action = tool_bar->addAction("n-", graph_area, SLOT(decrease_triangulation()));
    action->setShortcut(QString("5"));
    
    action = tool_bar->addAction("p++", graph_area, SLOT(increase_protrusion()));
    action->setShortcut(QString("6"));
    
    action = tool_bar->addAction("p--", graph_area, SLOT(decrease_protrusion()));
    action->setShortcut(QString("7"));
	
	action = tool_bar->addAction("m+", graph_area, SLOT(increase_m()));
    action->setShortcut(QString("8"));
    
    action = tool_bar->addAction("m-", graph_area, SLOT(decrease_m()));
    action->setShortcut(QString("9"));

    action = tool_bar->addAction("Exit", graph_area, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);

    // window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");
    window->show();
    app.exec();
    
    delete graph_area;
    delete tool_bar;
    delete window;
    return 0;
}
