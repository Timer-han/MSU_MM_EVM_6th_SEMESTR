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

    action = tool_bar->addAction("Change function", graph_area, SLOT(ChangeFunc()));
    action->setShortcut(QString("0"));
    
    action = tool_bar->addAction("Change Type Of Graph", graph_area, SLOT(ChangeTypeOfGraph()));
    action->setShortcut(QString("1"));
    
    action = tool_bar->addAction("Compress area", graph_area, SLOT(CompressArea()));
    action->setShortcut(QString("2"));
    
    action = tool_bar->addAction("Expand area", graph_area, SLOT(ExtendArea()));
    action->setShortcut(QString("3"));
    
    action = tool_bar->addAction("Increase N", graph_area, SLOT(IncreaseN()));
    action->setShortcut(QString("4"));
    
    action = tool_bar->addAction("Decrease N", graph_area, SLOT(DecreaseN()));
    action->setShortcut(QString("5"));
    
    action = tool_bar->addAction("Increase Func Middle", graph_area, SLOT(IncreaseFuncMiddle()));
    action->setShortcut(QString("6"));
    
    action = tool_bar->addAction("Decrease Func Middle", graph_area, SLOT(DecreaseFuncMiddle()));
    action->setShortcut(QString("7"));
	
	action = tool_bar->addAction("Increase M", graph_area, SLOT(IncreaseM()));
    action->setShortcut(QString("8"));
    
    action = tool_bar->addAction("Decrease M", graph_area, SLOT(DecreaseM()));
    action->setShortcut(QString("9"));

    action = tool_bar->addAction("Exit", graph_area, SLOT(Finish()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");
    window->show();
    app.exec();
    
    delete graph_area;
    delete tool_bar;
    delete window;
    return 0;
}
