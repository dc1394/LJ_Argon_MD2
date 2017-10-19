/*! \file Ar_moleculardynamics.h
\brief アルゴンに対して、分子動力学シミュレーションを行うクラスの実装

Copyright ©  2015 @dc1394 All Rights Reserved.
This software is released under the BSD 2-Clause License.
*/

#include "Ar_moleculardynamics.h"
#include "myrandom/myrand.h"
#include <cmath>                    // for std::sqrt, std::pow
#include <random>                   // for std::uniform_real_distribution
#include <boost/assert.hpp>         // for BOOST_ASSERT
#include <tbb/combinable.h>         // for tbb::combinable
#include <tbb/parallel_for.h>       // for tbb::parallel_for

namespace moleculardynamics {
    // #region static private 定数

    double const Ar_moleculardynamics::FIRSTSCALE = 1.0;

    double const Ar_moleculardynamics::FIRSTTEMP = 50.0;

    double const Ar_moleculardynamics::SIGMA = 3.405E-10;

    double const Ar_moleculardynamics::VDW_RADIUS = 1.88E-10;

    double const Ar_moleculardynamics::ALPHA = 0.2;

    double const Ar_moleculardynamics::ATM = 9.86923266716013E-6;

    double const Ar_moleculardynamics::AVOGADRO_CONSTANT = 6.022140857E+23;

    double const Ar_moleculardynamics::DT = 0.0001;

    double const Ar_moleculardynamics::HARTREE = 4.35974465054E-18;

    double const Ar_moleculardynamics::KB = 1.3806488E-23;

    double const Ar_moleculardynamics::TAU =
        std::sqrt(0.039948 / Ar_moleculardynamics::AVOGADRO_CONSTANT * Ar_moleculardynamics::SIGMA * Ar_moleculardynamics::SIGMA / Ar_moleculardynamics::YPSILON);

    double const Ar_moleculardynamics::YPSILON = 1.6540172624E-21;

    // #endregion static private 定数

    // #region コンストラクタ

    Ar_moleculardynamics::Ar_moleculardynamics()
        :
        Atoms([this] { return std::cref(atoms_); }, nullptr),
        MD_iter([this] { return MD_iter_; }, nullptr),
        Nc([this] { return Nc_; }, nullptr),
        NumAtom([this] { return NumAtom_; }, nullptr),
        periodiclen([this] { return periodiclen_; }, nullptr),
        Uk([this] { return DimensionlessToHartree(Uk_); }, nullptr),
        Up([this] { return DimensionlessToHartree(Up_); }, nullptr),
        Utot([this] { return DimensionlessToHartree(Utot_); }, nullptr),
        atoms_(Nc_ * Nc_ * Nc_ * 4),
        dt2(DT * DT),
        rc2_(SystemParam::RCUTOFF * SystemParam::RCUTOFF),
        rcm6_(std::pow(SystemParam::RCUTOFF, -6.0)),
        rcm12_(std::pow(SystemParam::RCUTOFF, -12.0)),
        Tg_(Ar_moleculardynamics::FIRSTTEMP * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON),
        Vrc_(4.0 * (rcm6_ - rcm12_))
    {
        // initalize parameters
        lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;

        recalc();
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    void Ar_moleculardynamics::runCalc()
    {
        moveAtoms();
        checkPairlist();
        calcForcePair();
        moveAtoms();
        periodic();

        // 繰り返し回数と時間を増加
        t_ = static_cast<double>(MD_iter_)* Ar_moleculardynamics::DT;
        MD_iter_++;
    }

    double Ar_moleculardynamics::getDeltat() const
    {
        return Ar_moleculardynamics::TAU * t_ * 1.0E+12;
    }

    float Ar_moleculardynamics::getForce(std::int32_t n) const
    {
        return static_cast<float>(atoms_[n].f.norm());
    }

    double Ar_moleculardynamics::getLatticeconst() const
    {
        return Ar_moleculardynamics::SIGMA * lat_ * 1.0E+9;
    }

    double Ar_moleculardynamics::getPeriodiclen() const
    {
        return Ar_moleculardynamics::SIGMA * periodiclen_ * 1.0E+9;
    }

    double Ar_moleculardynamics::getPressure() const
    {
        auto const V = std::pow(Ar_moleculardynamics::SIGMA * periodiclen_, 3);
        auto const ideal = NumAtom * Ar_moleculardynamics::YPSILON * Tc_;

        return (ideal - virial_ * Ar_moleculardynamics::YPSILON / 3.0) / V * Ar_moleculardynamics::ATM;
    }

    double Ar_moleculardynamics::getTcalc() const
    {
        return Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB * Tc_;
    }

    double Ar_moleculardynamics::getTgiven() const
    {
        return Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::KB * Tg_;
    }

    void Ar_moleculardynamics::recalc()
    {
        t_ = 0.0;
        MD_iter_ = 1;

        MD_initPos();

        MD_initVel();

        periodiclen_ = lat_ * static_cast<double>(Nc_);

        m_ = static_cast<std::int32_t>(periodiclen_ / (SystemParam::RCUTOFF + SystemParam::MARGIN)) - 1;
        
        if (m_ > 2) {
            pmesh_.reset(new MeshList(periodiclen_));
            pmesh_->set_number_of_atoms(atoms_.size());
            pmesh_->make_pair(atoms_, pairs_);
        }
        else {
            makePair();
        }
    }

    void Ar_moleculardynamics::setEnsemble(EnsembleType ensemble)
    {
        ensemble_ = ensemble;
        recalc();
    }

    void Ar_moleculardynamics::setNc(std::int32_t Nc)
    {
        Nc_ = Nc;
        atoms_.resize(Nc_ * Nc_ * Nc_ * 4);

        ModLattice();
    }

    void Ar_moleculardynamics::setScale(double scale)
    {
        scale_ = scale;
        ModLattice();
    }

    void Ar_moleculardynamics::setTgiven(double Tgiven)
    {
        Tg_ = Tgiven * Ar_moleculardynamics::KB / Ar_moleculardynamics::YPSILON;
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void Ar_moleculardynamics::calcForcePair()
    {
        // 各原子に働く力の初期化
        for (auto && a : atoms_) {
            a.f = Eigen::Vector4d::Zero();
        }

        tbb::combinable<double> Up;
        tbb::combinable<double> virial;

        tbb::parallel_for(
            tbb::blocked_range<std::int32_t>(0, pairs_.size()),
            [this, &Up, &virial](tbb::blocked_range<std::int32_t> const & range) {
            for (auto && k = range.begin(); k != range.end(); ++k) {
                auto const i = pairs_[k].first;
                auto const j = pairs_[k].second;
                Eigen::Vector4d d = atoms_[j].r - atoms_[i].r;

                SystemParam::adjust_periodic(d, periodiclen_);
                auto const r2 = d.squaredNorm();

                if (r2 <= rc2_) {
                    auto const r6 = r2 * r2 * r2;
                    auto const dFdr = (24.0 * r6 - 48.0) / (r6 * r6 * r2);

                    atoms_[i].f += dFdr * d;
                    atoms_[j].f -= dFdr * d;

                    auto const df = dFdr * DT;

                    atoms_[i].p += df * d;
                    atoms_[j].p -= df * d;

                    auto const r12 = r6 * r6;
                    Up.local() += 4.0 * (1.0 / r12 - 1.0 / r6) + Vrc_;
                    virial.local() += r2 * dFdr;
                }
            }
        });

        Up_ = Up.combine(std::plus<double>());
        virial_ = virial.combine(std::plus<double>());
    }

    void Ar_moleculardynamics::checkPairlist()
    {
        auto vmax2 = 0.0;

        for (auto && a : atoms_) {
            auto const v2 = a.p.squaredNorm();
            if (vmax2 < v2) {
                vmax2 = v2;
            }
        }

        auto const vmax = std::sqrt(vmax2);
        margin_length_ -= vmax * 2.0 * DT;

        if (margin_length_ < 0.0) {
            margin_length_ = SystemParam::MARGIN;
            if (m_ > 2) {
                pmesh_->make_pair(atoms_, pairs_);
            }
            else {
                makePair();
            }
        }
    }

    double Ar_moleculardynamics::DimensionlessToHartree(double e) const
    {
        return e * Ar_moleculardynamics::YPSILON / Ar_moleculardynamics::HARTREE;
    }

    void Ar_moleculardynamics::makePair()
    {
        pairs_.clear();

        tbb::parallel_for(
            tbb::blocked_range<std::int32_t>(0, NumAtom_ - 1),
            [this](tbb::blocked_range<std::int32_t> const & range) {
            for (auto && i = range.begin(); i != range.end(); ++i) {
                for (auto j = i + 1; j < NumAtom_; j++) {
                    Eigen::Vector4d d = atoms_[j].r - atoms_[i].r;

                    SystemParam::adjust_periodic(d, periodiclen_);

                    if (d.squaredNorm() <= rc2_) {
                        pairs_.push_back(std::make_pair(i, j));
                    }
                }
            }
        });
    }

    void Ar_moleculardynamics::MD_initPos()
    {
        double sx, sy, sz;
        auto n = 0;

        for (auto i = 0; i < Nc_; i++) {
            for (auto j = 0; j < Nc_; j++) {
                for (auto k = 0; k < Nc_; k++) {
                    // 基本セルをコピーする
                    sx = static_cast<double>(i)* lat_;
                    sy = static_cast<double>(j)* lat_;
                    sz = static_cast<double>(k)* lat_;

                    // 基本セル内には4つの原子がある
                    atoms_[n].r = Eigen::Vector4d(sx, sy, sz, 0.0);
                    n++;

                    atoms_[n].r = Eigen::Vector4d(0.5 * lat_ + sx, 0.5 * lat_ + sy, sz, 0.0);
                    n++;

                    atoms_[n].r = Eigen::Vector4d(sx, 0.5 * lat_ + sy, 0.5 * lat_ + sz, 0.0);
                    n++;

                    atoms_[n].r = Eigen::Vector4d(0.5 * lat_ + sx, sy, 0.5 * lat_ + sz, 0.0);
                    n++;
                }
            }
        }

        NumAtom_ = n;

        // move the center of mass to the origin
        // 系の重心を座標系の原点とする
        sx = 0.0;
        sy = 0.0;
        sz = 0.0;

        for (auto n = 0; n < NumAtom_; n++) {
            sx += atoms_[n].r[0];
            sy += atoms_[n].r[1];
            sz += atoms_[n].r[2];
        }

        sx /= static_cast<double>(NumAtom_);
        sy /= static_cast<double>(NumAtom_);
        sz /= static_cast<double>(NumAtom_);

        for (auto n = 0; n < NumAtom_; n++) {
            atoms_[n].r -= Eigen::Vector4d(sx, sy, sz, 0.0);
        }
    }

    void Ar_moleculardynamics::MD_initVel()
    {
        auto const v = std::sqrt(3.0 * Tg_);

        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        myrandom::MyRand<std::uniform_real_distribution<double> > mr(dist);

        for (auto && a : atoms_) {
            Eigen::Vector4d rnd(mr.myrand(), mr.myrand(), mr.myrand(), 0.0);

            // 方向はランダムに与える
            a.p = v * rnd / rnd.norm();
        }

        Eigen::Vector4d s(0.0, 0.0, 0.0, 0.0);

        for (auto && a : atoms_) {
            s += a.p;
        }

        s /= static_cast<double>(NumAtom_);

        // 重心の並進運動を避けるために、速度の和がゼロになるように補正
        for (auto && a : atoms_) {
            a.p -= s;
        }
    }

    void Ar_moleculardynamics::ModLattice()
    {
        lat_ = std::pow(2.0, 2.0 / 3.0) * scale_;
        recalc();
    }

    void Ar_moleculardynamics::moveAtoms()
    {
        // 運動エネルギーの初期化
        Uk_ = 0.0;

        // calculate temperture
        for (auto && a : atoms_) {
            Uk_ += a.p.squaredNorm();
        }

        // 運動エネルギーの計算
        Uk_ *= 0.5;

        // 全エネルギー（運動エネルギー+ポテンシャルエネルギー）の計算
        Utot_ = Uk_ + Up_;

        // 温度の計算
        Tc_ = Uk_ / (1.5 * static_cast<double>(NumAtom_));

        switch (ensemble_) {
        case EnsembleType::NVE:
            {
                tbb::parallel_for(
                    tbb::blocked_range<std::int32_t>(0, NumAtom_),
                    [this](tbb::blocked_range<std::int32_t> const & range) {
                    for (auto && n = range.begin(); n != range.end(); ++n) {
                        atoms_[n].r += atoms_[n].p * DT * 0.5;
                    }
                });
            }   
            break;

        case EnsembleType::NVT:
            {
                //auto const s = std::sqrt((Tg_ + Ar_moleculardynamics::ALPHA * (Tc_ - Tg_)) / Tc_);

                //tbb::parallel_for(
                //    tbb::blocked_range<std::int32_t>(0, NumAtom_),
                //    [this, s](tbb::blocked_range<std::int32_t> const & range) {
                //    for (auto && n = range.begin(); n != range.end(); ++n) {
                //        atoms_[n].p *= s;
                //        atoms_[n].r += atoms_[n].p * DT * 0.5;
                //    }
                //});
                
                auto const gamma = 1.0;
                auto const D = std::sqrt(2.0 * gamma * Tg_ / DT);

                std::normal_distribution<double> nd(0.0, D);
                myrandom::MyRand<std::normal_distribution<double> > mr(nd);
                for (auto && atom : atoms_) {
                    atom.p[0] += (-gamma * atom.p[0] + mr.myrand()) * DT;
                    atom.p[1] += (-gamma * atom.p[1] + mr.myrand()) * DT;
                    atom.p[2] += (-gamma * atom.p[2] + mr.myrand()) * DT;

                    atom.r += atom.p * DT * 0.5;
                }
            }
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
        }
    }

    void Ar_moleculardynamics::periodic()
    {
        // consider the periodic boundary condination
        // セルの外側に出たら座標をセル内に戻す
        tbb::parallel_for(
            tbb::blocked_range<std::int32_t>(0, NumAtom_),
            [this](tbb::blocked_range<std::int32_t> const & range) {
            for (auto && n = range.begin(); n != range.end(); ++n) {
                if (atoms_[n].r[0] > periodiclen_) {
                    atoms_[n].r[0] -= periodiclen_;
                }
                else if (atoms_[n].r[0] < 0.0) {
                    atoms_[n].r[0] += periodiclen_;
                }

                if (atoms_[n].r[1] > periodiclen_) {
                    atoms_[n].r[1] -= periodiclen_;
                }
                else if (atoms_[n].r[1] < 0.0) {
                    atoms_[n].r[1] += periodiclen_;
                }

                if (atoms_[n].r[2] > periodiclen_) {
                    atoms_[n].r[2] -= periodiclen_;
                }
                else if (atoms_[n].r[2] < 0.0) {
                    atoms_[n].r[2] += periodiclen_;
                }
            }
        });
    }

    // #endregion privateメンバ関数
}
