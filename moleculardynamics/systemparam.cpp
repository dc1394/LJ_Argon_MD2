/*! \file systemparam.cpp
    \brief メッシュリストクラスの実装

    Copyright © 2017 @dc1394 All Rights Reserved.
    (but this is originally adapted by @kaityo256 for meshlist.cpp from https://github.com/kaityo256/mdstep/tree/master/step3 )
    This software is released under the BSD 2-Clause License.
*/

#include "systemparam.h"

namespace moleculardynamics {
    // #region publicメンバ変数 

    double const SystemParam::MARGIN = 0.75;
    
    double const SystemParam::ML2 = (SystemParam::RCUTOFF + SystemParam::MARGIN) * (SystemParam::RCUTOFF + SystemParam::MARGIN);
    
    double const SystemParam::RCUTOFF = 2.5;

    // #endregion publicメンバ変数
}
