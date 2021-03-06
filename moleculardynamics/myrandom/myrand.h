﻿/*! \file myrand.h
    \brief 自作乱数クラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _MYRAND_H_
#define _MYRAND_H_

#pragma once

#include <cstdint>                      // for std::uint_least32_t
#include <random>                       // for std::mt19937
#include <vector>                       // for std::vector
#include <boost/range/algorithm.hpp>    // for boost::generate

namespace myrandom {
    template <typename Distribution>
    //! A template class.
    /*!
        自作乱数クラス
        \tparam Distribution 乱数分布クラス
    */
    class MyRand final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param distribution 乱数の分布
        */
        explicit MyRand(Distribution const & distribution);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~MyRand() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //!  A public member function.
        /*!
            [min, max]の閉区間で一様乱数を生成する
        */
        double myrand()
        {
            return distribution_(randengine_);
        }

        // #endregion メンバ関数

        // #region メンバ変数

    private:
        //! A private static member variable (constant expression).
        /*!
            初期乱数生成用のstd::vectorのサイズ
        */
        static std::vector<std::uint_least32_t>::size_type const SIZE = 8;

        //! A private member variable.
        /*!
            乱数の分布
        */
        Distribution distribution_;
        
        //! A private member variable.
        /*!
            乱数エンジン
        */
        std::mt19937 randengine_;

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
        デフォルトコンストラクタ（禁止）
        */
        MyRand() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        MyRand(const MyRand &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        MyRand & operator=(const MyRand &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    template <typename Distribution>
    MyRand<Distribution>::MyRand(Distribution const & distribution)
        : distribution_(distribution)
    {
        // ランダムデバイス
        std::random_device rnd;

        // 初期化用ベクタ
        std::vector<std::uint_least32_t> v(SIZE);

        // ベクタの初期化
        boost::generate(v, std::ref(rnd));

        std::seed_seq seq(v.begin(), v.end());

        // 乱数エンジン
        randengine_ = std::mt19937(seq);
    }
}

#endif  // _MYRAND_H_
