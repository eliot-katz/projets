/****************************************************************************
** Meta object code from reading C++ file 'McWorker.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../app/gui/src/McWorker.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'McWorker.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_gui__McWorker_t {
    QByteArrayData data[49];
    char stringdata0[441];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_gui__McWorker_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_gui__McWorker_t qt_meta_stringdata_gui__McWorker = {
    {
QT_MOC_LITERAL(0, 0, 13), // "gui::McWorker"
QT_MOC_LITERAL(1, 14, 8), // "progress"
QT_MOC_LITERAL(2, 23, 0), // ""
QT_MOC_LITERAL(3, 24, 11), // "std::size_t"
QT_MOC_LITERAL(4, 36, 5), // "nDone"
QT_MOC_LITERAL(5, 42, 4), // "mean"
QT_MOC_LITERAL(6, 47, 11), // "halfwidth95"
QT_MOC_LITERAL(7, 59, 8), // "finished"
QT_MOC_LITERAL(8, 68, 5), // "price"
QT_MOC_LITERAL(9, 74, 9), // "std_error"
QT_MOC_LITERAL(10, 84, 6), // "ci_low"
QT_MOC_LITERAL(11, 91, 7), // "ci_high"
QT_MOC_LITERAL(12, 99, 11), // "n_effective"
QT_MOC_LITERAL(13, 111, 10), // "elapsed_ms"
QT_MOC_LITERAL(14, 122, 6), // "failed"
QT_MOC_LITERAL(15, 129, 3), // "why"
QT_MOC_LITERAL(16, 133, 8), // "canceled"
QT_MOC_LITERAL(17, 142, 14), // "greeksFinished"
QT_MOC_LITERAL(18, 157, 5), // "delta"
QT_MOC_LITERAL(19, 163, 4), // "vega"
QT_MOC_LITERAL(20, 168, 5), // "gamma"
QT_MOC_LITERAL(21, 174, 3), // "rho"
QT_MOC_LITERAL(22, 178, 5), // "theta"
QT_MOC_LITERAL(23, 184, 14), // "stressFinished"
QT_MOC_LITERAL(24, 199, 9), // "basePrice"
QT_MOC_LITERAL(25, 209, 13), // "stressedPrice"
QT_MOC_LITERAL(26, 223, 3), // "pnl"
QT_MOC_LITERAL(27, 227, 10), // "runPricing"
QT_MOC_LITERAL(28, 238, 22), // "rw::market::MarketData"
QT_MOC_LITERAL(29, 261, 3), // "mkt"
QT_MOC_LITERAL(30, 265, 22), // "rw::market::Instrument"
QT_MOC_LITERAL(31, 288, 4), // "inst"
QT_MOC_LITERAL(32, 293, 20), // "rw::config::McConfig"
QT_MOC_LITERAL(33, 314, 3), // "cfg"
QT_MOC_LITERAL(34, 318, 5), // "sigma"
QT_MOC_LITERAL(35, 324, 9), // "runGreeks"
QT_MOC_LITERAL(36, 334, 2), // "mc"
QT_MOC_LITERAL(37, 337, 24), // "rw::config::GreeksConfig"
QT_MOC_LITERAL(38, 362, 2), // "gk"
QT_MOC_LITERAL(39, 365, 6), // "method"
QT_MOC_LITERAL(40, 372, 9), // "runStress"
QT_MOC_LITERAL(41, 382, 7), // "baseMkt"
QT_MOC_LITERAL(42, 390, 8), // "baseInst"
QT_MOC_LITERAL(43, 399, 9), // "baseSigma"
QT_MOC_LITERAL(44, 409, 6), // "dS_rel"
QT_MOC_LITERAL(45, 416, 6), // "dSigma"
QT_MOC_LITERAL(46, 423, 2), // "dR"
QT_MOC_LITERAL(47, 426, 2), // "dT"
QT_MOC_LITERAL(48, 429, 11) // "requestStop"

    },
    "gui::McWorker\0progress\0\0std::size_t\0"
    "nDone\0mean\0halfwidth95\0finished\0price\0"
    "std_error\0ci_low\0ci_high\0n_effective\0"
    "elapsed_ms\0failed\0why\0canceled\0"
    "greeksFinished\0delta\0vega\0gamma\0rho\0"
    "theta\0stressFinished\0basePrice\0"
    "stressedPrice\0pnl\0runPricing\0"
    "rw::market::MarketData\0mkt\0"
    "rw::market::Instrument\0inst\0"
    "rw::config::McConfig\0cfg\0sigma\0runGreeks\0"
    "mc\0rw::config::GreeksConfig\0gk\0method\0"
    "runStress\0baseMkt\0baseInst\0baseSigma\0"
    "dS_rel\0dSigma\0dR\0dT\0requestStop"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_gui__McWorker[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       6,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    3,   64,    2, 0x06 /* Public */,
       7,    6,   71,    2, 0x06 /* Public */,
      14,    1,   84,    2, 0x06 /* Public */,
      16,    0,   87,    2, 0x06 /* Public */,
      17,    6,   88,    2, 0x06 /* Public */,
      23,    4,  101,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      27,    4,  110,    2, 0x0a /* Public */,
      35,    6,  119,    2, 0x0a /* Public */,
      40,    8,  132,    2, 0x0a /* Public */,
      48,    0,  149,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3, QMetaType::Double, QMetaType::Double,    4,    5,    6,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double, 0x80000000 | 3, QMetaType::LongLong,    8,    9,   10,   11,   12,   13,
    QMetaType::Void, QMetaType::QString,   15,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::LongLong,   18,   19,   20,   21,   22,   13,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::LongLong,   24,   25,   26,   13,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 28, 0x80000000 | 30, 0x80000000 | 32, QMetaType::Double,   29,   31,   33,   34,
    QMetaType::Void, 0x80000000 | 28, 0x80000000 | 30, 0x80000000 | 32, 0x80000000 | 37, QMetaType::Double, QMetaType::QString,   29,   31,   36,   38,   34,   39,
    QMetaType::Void, 0x80000000 | 28, 0x80000000 | 30, 0x80000000 | 32, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double,   41,   42,   33,   43,   44,   45,   46,   47,
    QMetaType::Void,

       0        // eod
};

void gui::McWorker::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<McWorker *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->progress((*reinterpret_cast< std::size_t(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 1: _t->finished((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4])),(*reinterpret_cast< std::size_t(*)>(_a[5])),(*reinterpret_cast< long long(*)>(_a[6]))); break;
        case 2: _t->failed((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 3: _t->canceled(); break;
        case 4: _t->greeksFinished((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4])),(*reinterpret_cast< double(*)>(_a[5])),(*reinterpret_cast< long long(*)>(_a[6]))); break;
        case 5: _t->stressFinished((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< long long(*)>(_a[4]))); break;
        case 6: _t->runPricing((*reinterpret_cast< rw::market::MarketData(*)>(_a[1])),(*reinterpret_cast< rw::market::Instrument(*)>(_a[2])),(*reinterpret_cast< rw::config::McConfig(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4]))); break;
        case 7: _t->runGreeks((*reinterpret_cast< rw::market::MarketData(*)>(_a[1])),(*reinterpret_cast< rw::market::Instrument(*)>(_a[2])),(*reinterpret_cast< rw::config::McConfig(*)>(_a[3])),(*reinterpret_cast< rw::config::GreeksConfig(*)>(_a[4])),(*reinterpret_cast< double(*)>(_a[5])),(*reinterpret_cast< QString(*)>(_a[6]))); break;
        case 8: _t->runStress((*reinterpret_cast< rw::market::MarketData(*)>(_a[1])),(*reinterpret_cast< rw::market::Instrument(*)>(_a[2])),(*reinterpret_cast< rw::config::McConfig(*)>(_a[3])),(*reinterpret_cast< double(*)>(_a[4])),(*reinterpret_cast< double(*)>(_a[5])),(*reinterpret_cast< double(*)>(_a[6])),(*reinterpret_cast< double(*)>(_a[7])),(*reinterpret_cast< double(*)>(_a[8]))); break;
        case 9: _t->requestStop(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (McWorker::*)(std::size_t , double , double );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::progress)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (McWorker::*)(double , double , double , double , std::size_t , long long );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::finished)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (McWorker::*)(QString );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::failed)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (McWorker::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::canceled)) {
                *result = 3;
                return;
            }
        }
        {
            using _t = void (McWorker::*)(double , double , double , double , double , long long );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::greeksFinished)) {
                *result = 4;
                return;
            }
        }
        {
            using _t = void (McWorker::*)(double , double , double , long long );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&McWorker::stressFinished)) {
                *result = 5;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject gui::McWorker::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_gui__McWorker.data,
    qt_meta_data_gui__McWorker,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *gui::McWorker::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *gui::McWorker::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_gui__McWorker.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int gui::McWorker::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}

// SIGNAL 0
void gui::McWorker::progress(std::size_t _t1, double _t2, double _t3)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void gui::McWorker::finished(double _t1, double _t2, double _t3, double _t4, std::size_t _t5, long long _t6)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t4))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t5))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t6))) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void gui::McWorker::failed(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void gui::McWorker::canceled()
{
    QMetaObject::activate(this, &staticMetaObject, 3, nullptr);
}

// SIGNAL 4
void gui::McWorker::greeksFinished(double _t1, double _t2, double _t3, double _t4, double _t5, long long _t6)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t4))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t5))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t6))) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void gui::McWorker::stressFinished(double _t1, double _t2, double _t3, long long _t4)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t4))) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
