/****************************************************************************
** Meta object code from reading C++ file 'MainWindow.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../app/gui/src/MainWindow.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MainWindow.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[27];
    char stringdata0[243];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 7), // "onRunMC"
QT_MOC_LITERAL(2, 19, 0), // ""
QT_MOC_LITERAL(3, 20, 11), // "onRunGreeks"
QT_MOC_LITERAL(4, 32, 11), // "onRunStress"
QT_MOC_LITERAL(5, 44, 14), // "onLoadDefaults"
QT_MOC_LITERAL(6, 59, 7), // "onReset"
QT_MOC_LITERAL(7, 67, 12), // "onMcProgress"
QT_MOC_LITERAL(8, 80, 11), // "std::size_t"
QT_MOC_LITERAL(9, 92, 1), // "n"
QT_MOC_LITERAL(10, 94, 4), // "mean"
QT_MOC_LITERAL(11, 99, 4), // "half"
QT_MOC_LITERAL(12, 104, 12), // "onMcFinished"
QT_MOC_LITERAL(13, 117, 5), // "price"
QT_MOC_LITERAL(14, 123, 2), // "se"
QT_MOC_LITERAL(15, 126, 2), // "lo"
QT_MOC_LITERAL(16, 129, 2), // "hi"
QT_MOC_LITERAL(17, 132, 4), // "nEff"
QT_MOC_LITERAL(18, 137, 2), // "ms"
QT_MOC_LITERAL(19, 140, 10), // "onMcFailed"
QT_MOC_LITERAL(20, 151, 3), // "why"
QT_MOC_LITERAL(21, 155, 12), // "onMcCanceled"
QT_MOC_LITERAL(22, 168, 14), // "onSnapBaseline"
QT_MOC_LITERAL(23, 183, 13), // "onSaveProject"
QT_MOC_LITERAL(24, 197, 13), // "onLoadProject"
QT_MOC_LITERAL(25, 211, 14), // "onShortcutSave"
QT_MOC_LITERAL(26, 226, 16) // "onShortcutSaveAs"

    },
    "MainWindow\0onRunMC\0\0onRunGreeks\0"
    "onRunStress\0onLoadDefaults\0onReset\0"
    "onMcProgress\0std::size_t\0n\0mean\0half\0"
    "onMcFinished\0price\0se\0lo\0hi\0nEff\0ms\0"
    "onMcFailed\0why\0onMcCanceled\0onSnapBaseline\0"
    "onSaveProject\0onLoadProject\0onShortcutSave\0"
    "onShortcutSaveAs"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

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
       1,    0,   84,    2, 0x08 /* Private */,
       3,    0,   85,    2, 0x08 /* Private */,
       4,    0,   86,    2, 0x08 /* Private */,
       5,    0,   87,    2, 0x08 /* Private */,
       6,    0,   88,    2, 0x08 /* Private */,
       7,    3,   89,    2, 0x08 /* Private */,
      12,    6,   96,    2, 0x08 /* Private */,
      19,    1,  109,    2, 0x08 /* Private */,
      21,    0,  112,    2, 0x08 /* Private */,
      22,    0,  113,    2, 0x08 /* Private */,
      23,    0,  114,    2, 0x08 /* Private */,
      24,    0,  115,    2, 0x08 /* Private */,
      25,    0,  116,    2, 0x08 /* Private */,
      26,    0,  117,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 8, QMetaType::Double, QMetaType::Double,    9,   10,   11,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double, 0x80000000 | 8, QMetaType::LongLong,   13,   14,   15,   16,   17,   18,
    QMetaType::Void, QMetaType::QString,   20,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MainWindow *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->onRunMC(); break;
        case 1: _t->onRunGreeks(); break;
        case 2: _t->onRunStress(); break;
        case 3: _t->onLoadDefaults(); break;
        case 4: _t->onReset(); break;
        case 5: _t->onMcProgress((*reinterpret_cast< std::size_t(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 6: _t->onMcFinished((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4])),(*reinterpret_cast< std::size_t(*)>(_a[5])),(*reinterpret_cast< long long(*)>(_a[6]))); break;
        case 7: _t->onMcFailed((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: _t->onMcCanceled(); break;
        case 9: _t->onSnapBaseline(); break;
        case 10: _t->onSaveProject(); break;
        case 11: _t->onLoadProject(); break;
        case 12: _t->onShortcutSave(); break;
        case 13: _t->onShortcutSaveAs(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject MainWindow::staticMetaObject = { {
    QMetaObject::SuperData::link<QMainWindow::staticMetaObject>(),
    qt_meta_stringdata_MainWindow.data,
    qt_meta_data_MainWindow,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
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
