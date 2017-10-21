/*! \file systemparam.h
    \brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "systemparam.h"

namespace moleculardynamics {
    // #region static public 定数

    double const SystemParam::MARGIN = 0.75;
    
    double const SystemParam::ML2 = (SystemParam::RCUTOFF + SystemParam::MARGIN) * (SystemParam::RCUTOFF + SystemParam::MARGIN);
    
    double const SystemParam::RCUTOFF = 2.5;

    // #endregion static public 定数
}
