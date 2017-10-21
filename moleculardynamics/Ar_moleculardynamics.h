﻿/*! \file Ar_moleculardynamics.h
    \brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _AR_MOLECULARDYNAMICS_H_
#define _AR_MOLECULARDYNAMICS_H_

#pragma once

#include "../utility/property.h"
#include "meshlist.h"
#include "systemparam.h"
#include <cstdint>                  // for std::int32_t
#include <memory>                   // for std::unique_ptr

namespace moleculardynamics {
    using namespace utility;

    //! A enum.
    /*!
        アンサンブルのタイプの列挙型
    */
    enum class EnsembleType : std::int32_t {
        // NVEアンサンブル
        NVE = 0,

        // NVTアンサンブル
        NVT = 1
    };

    //! A enum.
    /*!
        温度制御の方法の列挙型
    */
    enum class TempControlType : std::int32_t {
        // Langevin法
        LANGEVIN = 0,

        // Nose-Hoover法
        NOSE_HOOVER = 1,

        // Woodcockの速度スケーリング法
        VELOCITY = 2
    };;

    //! A class.
    /*!
        アルゴンに対して、分子動力学シミュレーションを行うクラス
    */
    class Ar_moleculardynamics final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            コンストラクタ
        */
        Ar_moleculardynamics();

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Ar_moleculardynamics() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (constant).
        /*!
            シミュレーションを開始してからの経過時間を求める
        */
        double getDeltat() const;

        //! A public member function (constant).
        /*!
            n番目の原子に働く力を求める
        */
        float getForce(std::int32_t n) const;

        //! A public member function (constant).
        /*!
            格子定数を求める
        */
        double getLatticeconst() const;

        //! A public member function (constant).
        /*!
            周期境界条件の長さを求める
        */
        double getPeriodiclen() const;

        //! A public member function (constant).
        /*!
            計算された圧力を求める
        */
        double getPressure() const;

        //! A public member function (constant).
        /*!
            計算された温度の絶対温度を求める
        */
        double getTcalc() const;

        //! A public member function (constant).
        /*!
            与えた温度の絶対温度を求める
        */
        double getTgiven() const;

        //! A oublic member function.
        /*!
            再計算する
        */
        void recalc();

        //! A oublic member function.
        /*!
            MDを1ステップ計算する
        */
        void runCalc();

        //! A public member function.
        /*!
            アンサンブルを設定する
            \param ensemble 設定するアンサンブル
        */
        void setEnsemble(EnsembleType ensemble);

        //! A public member function.
        /*!
            スーパーセルの大きさを設定する
            \param Nc スーパーセルの大きさ
        */
        void setNc(std::int32_t Nc);

        //! A public member function.
        /*!
            格子定数のスケールを設定する
            \param scale 設定する格子定数のスケール
        */
        void setScale(double scale);

        //! A public member function.
        /*!
            温度を設定する
            \param Tgiven 設定する温度（絶対温度）
        */
        void setTgiven(double Tgiven);

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            原子に働く力を計算する
        */
        void calcForcePair();

        //! A private member function.
        /*!
            ペアリストの寿命をチェックする
        */
        void checkPairlist();

        //! A private member function.
        /*!
            エネルギーの単位を無次元単位からHartreeに変換する
            \param e 無次元単位で表されたエネルギー
            \return Hartree単位で表されたエネルギー
        */
        double DimensionlessToHartree(double e) const;

        //! A private member function.
        /*!
            Langevin法
        */
        void Langevin();

        //! A private member function.
        /*!
            ペアリストを構築する
        */
        void makePair();

        //! A private member function.
        /*!
            原子の初期位置を決める
        */
        void MD_initPos();

        //! A private member function.
        /*!
            原子の初期速度を決める
        */
        void MD_initVel();

        //! A private member function.
        /*!
            格子定数が変更されたときに呼ばれる
        */
        void ModLattice();

        //! A privte member function.
        /*!
            原子を移動させる
        */
        void moveAtoms();

        //! A private member function.
        /*!
            周期境界条件を用いて、原子の位置を補正する
        */
        void periodic();
        
        //! A private member function.
        /*!
            Woodcockの速度スケーリング法
        */
        void Woodcock_velocity_scaling();

        // #endregion privateメンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            原子へのプロパティ
        */
        Property<SystemParam::myatomvector const &> const Atoms;

        //! A property.
        /*!
            MDのステップ数へのプロパティ
        */
        Property<std::int32_t> const MD_iter;

        //! A property.
        /*!
            スーパーセルの個数へのプロパティ
        */
        Property<std::int32_t> const Nc;

        //! A property.
        /*!
            原子数へのプロパティ
        */
        Property<std::int32_t> const NumAtom;

        //! A property.
        /*!
            格子定数へのプロパティ
        */
        Property<double> const periodiclen;

        //! A property.
        /*!
            温度制御の方法へのプロパティ
        */
        Property<TempControlType> tempcon;

        //! A property.
        /*!
            運動エネルギーへのプロパティ
        */
        Property<double> const Uk;

        //! A property.
        /*!
            ポテンシャルエネルギーへのプロパティ
        */
        Property<double> const Up;

        //! A property.
        /*!
            全エネルギーへのプロパティ
        */
        Property<double> const Utot;

        // #endregion プロパティ

        // #region メンバ変数

    public:
        //! A public member variable (static constant).
        /*!
            初期のスーパーセルの個数
        */
        static auto const FIRSTNC = 6;

        //! A public member variable (static consttant).
        /*!
            初期の格子定数のスケール
        */
        static double const FIRSTSCALE;

        //! A public member variable (static constant).
        /*!
            初期温度（絶対温度）
        */
        static double const FIRSTTEMP;

        //! A public member variable (static constant).
        /*!
            アルゴン原子に対するσ
        */
        static double const SIGMA;

        //! A public member variable (static constant).
        /*!
            アルゴン原子のVan der Waals半径
        */
        static double const VDW_RADIUS;

    private:
        //! A private member variable (static constant).
        /*!
            Woodcockの温度スケーリングの係数
        */
        static double const ALPHA;

        //! A private member variable (static constant).
        /*!
            標準気圧
        */
        static double const ATM;

        //! A private member variable (static constant).
        /*!
            アボガドロ定数
        */
        static double const AVOGADRO_CONSTANT;

        //! A private member variable (static constant).
        /*!
            時間刻みΔt
        */
        static double const DT;

        //! A private member variable (static constant).
        /*!
            Langevin法の定数
        */
        static double const GAMMA;

        //! A private member variable (static constant).
        /*!
            1Hartree
        */
        static double const HARTREE;

        //! A private member variable (static constant).
        /*!
            ボルツマン定数
        */
        static double const KB;

        //! A private member variable (static constant).
        /*!
            アルゴン原子に対するτ
        */
        static double const TAU;

        //! A private member variable (static constant).
        /*!
            アルゴン原子に対するε
        */
        static double const YPSILON;

        //! A private member variable (static constant).
        /*!
            スーパーセルの個数
        */
        std::int32_t Nc_ = Ar_moleculardynamics::FIRSTNC;

        //! A private member variable.
        /*!
            原子の可変長配列
        */
        SystemParam::myatomvector atoms_;

        //! A private member variable (constant).
        /*!
            時間刻みの二乗
        */
        double const dt2;

        //! A private member variable.
        /*!
            アンサンブル
        */
        EnsembleType ensemble_ = EnsembleType::NVT;

        //! A private member variable.
        /*!
            格子定数
        */
        double lat_;
        
        //! A private member variable.
        /*!
            一辺のメッシュの数
        */
        std::int32_t m_;

        //! A private member variable.
        /*!
            ペアリストの寿命の長さ
        */
        double margin_length_;

        //! A private member variable.
        /*!
            MDのステップ数
        */
        std::int32_t MD_iter_;

        //! A private member variable (constant).
        /*!
            相互作用を計算するセルの個数
        */
        std::int32_t const ncp_ = 3;

        //! A private member variable.
        /*!
            メッシュのリストへのスマートポインタ
        */
        std::unique_ptr<MeshList> pmesh_;
        
        //! A private member variable.
        /*!
            原子数
        */
        std::int32_t NumAtom_;

        //! A private member variable.
        /*!
            ペアリスト
        */
        SystemParam::mypairvector pairs_;
        
        //! A private member variable.
        /*!
            周期境界条件の長さ
        */
        double periodiclen_;

        //! A private member variable (constant).
        /*!
            カットオフ半径の2乗
        */
        double const rc2_;

        //! A private member variable (constant).
        /*!
            カットオフ半径の逆数の6乗
        */
        double const rcm6_;

        //! A private member variable (constant).
        /*!
            カットオフ半径の逆数の12乗
        */
        double const rcm12_;

        //! A private member variable.
        /*!
            格子定数のスケーリングの定数
        */
        double scale_ = Ar_moleculardynamics::FIRSTSCALE;

        //! A private member variable.
        /*!
            時間
        */
        double t_;

        //! A private member variable.
        /*!
            計算された温度Tcalc
        */
        double Tc_;

        //! A private member variable.
        /*!
            温度制御の方法
        */
        TempControlType tempcon_ = TempControlType::VELOCITY;
        
        //! A private member variable.
        /*!
            与える温度Tgiven
        */
        double Tg_;

        //! A private member variable (constant).
        /*!
            運動エネルギー
        */
        double Uk_;

        //! A private member variable (constant).
        /*!
            ポテンシャルエネルギー
        */
        double Up_;

        //! A private member variable (constant).
        /*!
            全エネルギー
        */
        double Utot_;

        //! A private member variable (constant).
        /*!
            ビリアル
        */
        double virial_;

        //! A private member variable (constant).
        /*!
            ポテンシャルエネルギーの打ち切り
        */
        double const Vrc_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Ar_moleculardynamics(Ar_moleculardynamics const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Ar_moleculardynamics & operator=(Ar_moleculardynamics const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif      // _AR_MOLECULARDYNAMICS_H_
