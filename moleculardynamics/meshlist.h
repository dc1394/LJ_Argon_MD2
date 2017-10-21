/*! \file meshlist.h
    \brief メッシュリストクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    (but this is originally adapted by @kaityo256 for meshlist.hpp from https://github.com/kaityo256/mdstep/tree/master/step3 )
    This software is released under the BSD 2-Clause License.
*/

#ifndef _MESHLIST_H_
#define _MESHLIST_H_

#pragma once

#include "systemparam.h"

namespace moleculardynamics {
    //! A class.
    /*!
        メッシュリストクラス    
    */
    class MeshList final {
        // #region コンストラクタ・デストラクタ
        
    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param periodiclen 周期の長さ
        */
        explicit MeshList(double periodiclen);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~MeshList() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数
        
        //! A public member function.
        /*!
            原子の住所録を作成する
            \param atoms 原子の座標が格納された可変長配列
            \param pairs 原子のペアが格納された可変長配列
        */
        void make_pair(SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);
        
        //! A public member function.
        /*!
            原子の数を設定する
            \param pn 原子の数
        */
        void set_number_of_atoms(std::size_t pn) { sorted_buffer.resize(pn); }
        
        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            住所録から逆引きして調べる関数
            \param id 番地
            \param atoms 原子の座標が格納された可変長配列
            \param pairs 原子のペアが格納された可変長配列
        */
        void search(std::int32_t id, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);
        
        //! A private member function.
        /*!
            隣接番地での探索を行う関数
            \param id 番地
            \param ix 番地（x座標）
            \param iy 番地（y座標）
            \param iz 番地（z座標）
            \param atoms 原子の座標が格納された可変長配列
            \param pairs 原子のペアが格納された可変長配列
        */
        void search_other(std::int32_t id, std::int32_t ix, std::int32_t iy, std::int32_t iz, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);

        // #endregion privateメンバ関数

        // #region privateメンバ変数

        //! A private member variable.
        /*!
            どの番地に何個原子がいるかの数
        */
        std::vector<std::int32_t> count_;

        //! A private member variable.
        /*!
            番地番号でソートした際に、番地番号の頭出しのインデックス
        */
        std::vector<std::int32_t> indexes_;
        
        //! A private member variable.
        /*!
            一辺のメッシュの数    
        */
        std::int32_t m_;

        //! A private member variable.
        /*!
            メッシュのサイズ
        */
        double mesh_size_;

        //! A private member variable.
        /*!
            トータルのメッシュの数
        */
        std::int32_t number_of_mesh_;

        //! A private member variable.
        /*!
            番地番号でソートした原子インデックス
        */
        std::vector<std::int32_t> sorted_buffer;

        //! A private member variable (constant).
        /*!
            周期の長さ
        */
        double const periodiclen_;
        
        // #endregion privateメンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private copy constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        MeshList() = delete;
        
        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        MeshList(MeshList const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        MeshList & operator=(MeshList const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif	// _MESHLIST_H_
