/*! \file systemparam.h
    \brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _SYSTEMPARAM_H_
#define _SYSTEMPARAM_H_

#pragma once

#include <cstdint>                              // for std::int32_t
#include <utility>                              // for std::pair
#include <vector>                               // for std::vector
#include <Eigen/Core>                           // for Eigen::Vector4d
#include <boost/align/aligned_allocator.hpp>    // for boost::alignment::aligned_allocator

namespace moleculardynamics {
    //! A struct.
    /*!
        原子の情報が格納された構造体
    */
    #pragma pack(16)
    struct Atom {
        Eigen::Vector4d f;
        Eigen::Vector4d p;
        Eigen::Vector4d r;
    };

    //! A struct.
    /*!
        型エイリアスや定数が格納された構造体
    */
	struct SystemParam {
        // #region 型エイリアス

        using myatomvector = std::vector<Atom, boost::alignment::aligned_allocator<Atom> >;

        using mypairvector = std::vector<std::pair<std::int32_t, std::int32_t> >;

        // #endregion 型エイリアス

        // #region static publicメンバ関数

        //! A public static member function.
        /*!
            周期的境界条件の補正をする
            \param d x方向の補正
            \param periodiclen 周期の長さ
        */
        inline static void adjust_periodic(Eigen::Vector4d & d, double periodiclen);

        // #endregion static publicメンバ関数

        // #region publicメンバ変数

        //! A public member variable (static constant).
        /*!
            マージン
        */
        static double const MARGIN;

        //! A public member variable (static constant).
        /*!
            マージンと、マージンのカットオフの和の平方
        */
        static double const ML2;

        //! A public member variable (static constant).
        /*!
            カットオフ半径
        */
        static double const RCUTOFF;

        // #endregion publicメンバ変数
	};

    // #region publicメンバ関数の実装

    void SystemParam::adjust_periodic(Eigen::Vector4d & d, double periodiclen)
    {
        auto const LH = periodiclen * 0.5;

        if (d[0] < -LH) {
            d[0] += periodiclen;
        }
        else if (d[0] > LH) {
            d[0] -= periodiclen;
        }

        if (d[1] < -LH) {
            d[1] += periodiclen;
        }
        else if (d[1] > LH) {
            d[1] -= periodiclen;
        }

        if (d[2] < -LH) {
            d[2] += periodiclen;
        }
        else if (d[2] > LH) {
            d[2] -= periodiclen;
        }
    }

    // #endregion publicメンバ関数の実装
}

#endif	// _SYSTEMPARAM_H_