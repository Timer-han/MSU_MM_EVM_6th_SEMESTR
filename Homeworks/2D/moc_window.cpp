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
    QByteArrayData data[24];
    char stringdata0[340];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Window_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Window_t qt_meta_stringdata_Window = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Window"
QT_MOC_LITERAL(1, 7, 19), // "update_image_bounds"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 11), // "select_func"
QT_MOC_LITERAL(4, 40, 7), // "func_id"
QT_MOC_LITERAL(5, 48, 9), // "read_args"
QT_MOC_LITERAL(6, 58, 7), // "char*[]"
QT_MOC_LITERAL(7, 66, 4), // "argv"
QT_MOC_LITERAL(8, 71, 16), // "is_threads_ready"
QT_MOC_LITERAL(9, 88, 15), // "waiting_threads"
QT_MOC_LITERAL(10, 104, 14), // "memory_realloc"
QT_MOC_LITERAL(11, 119, 18), // "update_thread_data"
QT_MOC_LITERAL(12, 138, 17), // "save_prev_results"
QT_MOC_LITERAL(13, 156, 8), // "change_f"
QT_MOC_LITERAL(14, 165, 16), // "change_show_mode"
QT_MOC_LITERAL(15, 182, 21), // "increase_visible_area"
QT_MOC_LITERAL(16, 204, 21), // "decrease_visible_area"
QT_MOC_LITERAL(17, 226, 22), // "increase_triangulation"
QT_MOC_LITERAL(18, 249, 22), // "decrease_triangulation"
QT_MOC_LITERAL(19, 272, 19), // "increase_protrusion"
QT_MOC_LITERAL(20, 292, 19), // "decrease_protrusion"
QT_MOC_LITERAL(21, 312, 10), // "increase_m"
QT_MOC_LITERAL(22, 323, 10), // "decrease_m"
QT_MOC_LITERAL(23, 334, 5) // "close"

    },
    "Window\0update_image_bounds\0\0select_func\0"
    "func_id\0read_args\0char*[]\0argv\0"
    "is_threads_ready\0waiting_threads\0"
    "memory_realloc\0update_thread_data\0"
    "save_prev_results\0change_f\0change_show_mode\0"
    "increase_visible_area\0decrease_visible_area\0"
    "increase_triangulation\0decrease_triangulation\0"
    "increase_protrusion\0decrease_protrusion\0"
    "increase_m\0decrease_m\0close"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Window[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      19,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  109,    2, 0x0a /* Public */,
       3,    1,  110,    2, 0x0a /* Public */,
       5,    1,  113,    2, 0x0a /* Public */,
       8,    0,  116,    2, 0x0a /* Public */,
       9,    0,  117,    2, 0x0a /* Public */,
      10,    0,  118,    2, 0x0a /* Public */,
      11,    0,  119,    2, 0x0a /* Public */,
      12,    0,  120,    2, 0x0a /* Public */,
      13,    0,  121,    2, 0x0a /* Public */,
      14,    0,  122,    2, 0x0a /* Public */,
      15,    0,  123,    2, 0x0a /* Public */,
      16,    0,  124,    2, 0x0a /* Public */,
      17,    0,  125,    2, 0x0a /* Public */,
      18,    0,  126,    2, 0x0a /* Public */,
      19,    0,  127,    2, 0x0a /* Public */,
      20,    0,  128,    2, 0x0a /* Public */,
      21,    0,  129,    2, 0x0a /* Public */,
      22,    0,  130,    2, 0x0a /* Public */,
      23,    0,  131,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    4,
    QMetaType::Int, 0x80000000 | 6,    7,
    QMetaType::Bool,
    QMetaType::Void,
    QMetaType::Int,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void Window::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Window *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->update_image_bounds(); break;
        case 1: _t->select_func((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: { int _r = _t->read_args((*reinterpret_cast< char*(*)[]>(_a[1])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 3: { bool _r = _t->is_threads_ready();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 4: _t->waiting_threads(); break;
        case 5: { int _r = _t->memory_realloc();
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 6: _t->update_thread_data(); break;
        case 7: _t->save_prev_results(); break;
        case 8: _t->change_f(); break;
        case 9: _t->change_show_mode(); break;
        case 10: _t->increase_visible_area(); break;
        case 11: _t->decrease_visible_area(); break;
        case 12: _t->increase_triangulation(); break;
        case 13: _t->decrease_triangulation(); break;
        case 14: _t->increase_protrusion(); break;
        case 15: _t->decrease_protrusion(); break;
        case 16: _t->increase_m(); break;
        case 17: _t->decrease_m(); break;
        case 18: _t->close(); break;
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
        if (_id < 19)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 19;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 19)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 19;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
