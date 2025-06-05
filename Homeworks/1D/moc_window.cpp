/****************************************************************************
** Meta object code from reading C++ file 'window.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.13)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "window.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'window.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.13. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Window_t {
    QByteArrayData data[18];
    char stringdata0[196];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Window_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Window_t qt_meta_stringdata_Window = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Window"
QT_MOC_LITERAL(1, 7, 8), // "expand_n"
QT_MOC_LITERAL(2, 16, 0), // ""
QT_MOC_LITERAL(3, 17, 10), // "compress_n"
QT_MOC_LITERAL(4, 28, 10), // "increase_p"
QT_MOC_LITERAL(5, 39, 10), // "decrease_p"
QT_MOC_LITERAL(6, 50, 13), // "expand_bounds"
QT_MOC_LITERAL(7, 64, 15), // "compress_bounds"
QT_MOC_LITERAL(8, 80, 11), // "change_func"
QT_MOC_LITERAL(9, 92, 17), // "change_graph_type"
QT_MOC_LITERAL(10, 110, 13), // "paintLagrange"
QT_MOC_LITERAL(11, 124, 9), // "QPainter&"
QT_MOC_LITERAL(12, 134, 7), // "painter"
QT_MOC_LITERAL(13, 142, 3), // "pen"
QT_MOC_LITERAL(14, 146, 10), // "paintQubic"
QT_MOC_LITERAL(15, 157, 17), // "paintLagrResidual"
QT_MOC_LITERAL(16, 175, 2), // "pr"
QT_MOC_LITERAL(17, 178, 17) // "paintCubiResidual"

    },
    "Window\0expand_n\0\0compress_n\0increase_p\0"
    "decrease_p\0expand_bounds\0compress_bounds\0"
    "change_func\0change_graph_type\0"
    "paintLagrange\0QPainter&\0painter\0pen\0"
    "paintQubic\0paintLagrResidual\0pr\0"
    "paintCubiResidual"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Window[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   84,    2, 0x0a /* Public */,
       3,    0,   85,    2, 0x0a /* Public */,
       4,    0,   86,    2, 0x0a /* Public */,
       5,    0,   87,    2, 0x0a /* Public */,
       6,    0,   88,    2, 0x0a /* Public */,
       7,    0,   89,    2, 0x0a /* Public */,
       8,    0,   90,    2, 0x0a /* Public */,
       9,    0,   91,    2, 0x0a /* Public */,
      10,    2,   92,    2, 0x0a /* Public */,
      14,    2,   97,    2, 0x0a /* Public */,
      15,    3,  102,    2, 0x0a /* Public */,
      15,    2,  109,    2, 0x2a /* Public | MethodCloned */,
      17,    3,  114,    2, 0x0a /* Public */,
      17,    2,  121,    2, 0x2a /* Public | MethodCloned */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen,   12,   13,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen,   12,   13,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen, QMetaType::Bool,   12,   13,   16,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen,   12,   13,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen, QMetaType::Bool,   12,   13,   16,
    QMetaType::Void, 0x80000000 | 11, QMetaType::QPen,   12,   13,

       0        // eod
};

void Window::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Window *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->expand_n(); break;
        case 1: _t->compress_n(); break;
        case 2: _t->increase_p(); break;
        case 3: _t->decrease_p(); break;
        case 4: _t->expand_bounds(); break;
        case 5: _t->compress_bounds(); break;
        case 6: _t->change_func(); break;
        case 7: _t->change_graph_type(); break;
        case 8: _t->paintLagrange((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2]))); break;
        case 9: _t->paintQubic((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2]))); break;
        case 10: _t->paintLagrResidual((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2])),(*reinterpret_cast< bool(*)>(_a[3]))); break;
        case 11: _t->paintLagrResidual((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2]))); break;
        case 12: _t->paintCubiResidual((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2])),(*reinterpret_cast< bool(*)>(_a[3]))); break;
        case 13: _t->paintCubiResidual((*reinterpret_cast< QPainter(*)>(_a[1])),(*reinterpret_cast< QPen(*)>(_a[2]))); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Window::staticMetaObject = { {
    QMetaObject::SuperData::link<QWidget::staticMetaObject>(),
    qt_meta_stringdata_Window.data,
    qt_meta_data_Window,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Window::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Window::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Window.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int Window::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 14)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 14;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
