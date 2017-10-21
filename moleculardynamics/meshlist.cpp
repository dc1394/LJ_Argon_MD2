/*! \file meshlist.cpp
    \brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "meshlist.h"
#include <algorithm>        // for std::fill
#include <boost/assert.hpp> // for BOOST_ASSERT

namespace moleculardynamics {
    MeshList::MeshList(double periodiclen) : periodiclen_(periodiclen)
    {
        auto const SL = SystemParam::RCUTOFF + SystemParam::MARGIN;
        
        m_ = static_cast<std::int32_t>(periodiclen / SL) - 1;
        mesh_size_ = static_cast<double>(periodiclen) / m_;
        
        BOOST_ASSERT(m_ > 2);
        BOOST_ASSERT(mesh_size_ > SL);
        
        number_of_mesh_ = m_ * m_ * m_;
        count_.resize(number_of_mesh_);
        indexes_.resize(number_of_mesh_);
    }

    void MeshList::make_pair(SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs)
    {
        pairs.clear();
        
        auto const pn = atoms.size();

        std::vector<std::int32_t> particle_position(pn, 0);
        std::vector<std::int32_t> pointer(number_of_mesh_, 0);

        std::fill(count_.begin(), count_.end(), 0);

        auto const im = 1.0 / mesh_size_;
        for (auto i = 0; i < pn; i++) {
            auto ix = static_cast<std::int32_t>(atoms[i].r[0] * im);
            auto iy = static_cast<std::int32_t>(atoms[i].r[1] * im);
            auto iz = static_cast<std::int32_t>(atoms[i].r[2] * im);
            
            if (ix < 0) {
                ix += m_;
            }
            else if (ix >= m_) {
                ix -= m_;
            }

            if (iy < 0) {
                iy += m_;
            }
            else if (iy >= m_) {
                iy -= m_;
            }
            if (iz < 0) {
                iz += m_;
            }
            else if (iz >= m_) {
                iz -= m_;
            }

            auto const index = ix + iy * m_ + iz * m_ * m_;
            
            BOOST_ASSERT(index >= 0);
            BOOST_ASSERT(index < number_of_mesh_);
            
            count_[index]++;
            particle_position[i] = index;
        }
        
        indexes_[0] = 0;
        auto sum = 0;
        
        for (auto i = 0; i < number_of_mesh_ - 1; i++) {
            sum += count_[i];
            indexes_[i + 1] = sum;
        }
        
        for (auto i = 0; i < pn; i++) {
            auto const pos = particle_position[i];
            auto const j = indexes_[pos] + pointer[pos];
            sorted_buffer[j] = i;
            ++pointer[pos];
        }
        
        for (auto i = 0; i < number_of_mesh_; i++) {
            search(i, atoms, pairs);
        }
    }

    void MeshList::search_other(std::int32_t id, std::int32_t ix, std::int32_t iy, std::int32_t iz, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs)
    {
        if (ix < 0) {
            ix += m_;
        }
        else if (ix >= m_) {
            ix -= m_;
        }
        
        if (iy < 0) {
            iy += m_;
        }
        else if (iy >= m_) {
            iy -= m_;
        }
        
        if (iz < 0) {
            iz += m_;
        }
        else if (iz >= m_) {
            iz -= m_;
        }
        
        auto const id2 = ix + iy * m_ + iz * m_ * m_;

        for (auto k = indexes_[id]; k < indexes_[id] + count_[id]; k++) {
            for (auto m_ = indexes_[id2]; m_ < indexes_[id2] + count_[id2]; m_++) {
                auto const i = sorted_buffer[k];
                auto const j = sorted_buffer[m_];

                Eigen::Vector4d d = atoms[j].r - atoms[i].r;
                
                SystemParam::adjust_periodic(d, periodiclen_);
                
                if (d.squaredNorm() <= SystemParam::ML2) {
                    pairs.push_back(std::make_pair(i, j));
                }
            }
        }
    }

    void MeshList::search(std::int32_t id, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs)
    {
        auto const ix = id % m_;
        auto const iy = (id / m_) % m_;
        auto const iz = (id / m_ / m_);

        search_other(id, ix + 1, iy, iz, atoms, pairs);
        search_other(id, ix - 1, iy + 1, iz, atoms, pairs);
        search_other(id, ix, iy + 1, iz, atoms, pairs);
        search_other(id, ix + 1, iy + 1, iz, atoms, pairs);

        search_other(id, ix - 1, iy, iz + 1, atoms, pairs);
        search_other(id, ix, iy, iz + 1, atoms, pairs);
        search_other(id, ix + 1, iy, iz + 1, atoms, pairs);

        search_other(id, ix - 1, iy - 1, iz + 1, atoms, pairs);
        search_other(id, ix, iy - 1, iz + 1, atoms, pairs);
        search_other(id, ix + 1, iy - 1, iz + 1, atoms, pairs);

        search_other(id, ix - 1, iy + 1, iz + 1, atoms, pairs);
        search_other(id, ix, iy + 1, iz + 1, atoms, pairs);
        search_other(id, ix + 1, iy + 1, iz + 1, atoms, pairs);

        // Registration of self box
        auto const si = indexes_[id];
        auto const n = count_[id];

        for (auto k = si; k < si + n - 1; k++) {
            for (auto m_ = k + 1; m_ < si + n; m_++) {
                auto const i = sorted_buffer[k];
                auto const j = sorted_buffer[m_];

                Eigen::Vector4d d = atoms[j].r - atoms[i].r;

                SystemParam::adjust_periodic(d, periodiclen_);

                if (d.squaredNorm() <= SystemParam::ML2) {
                    pairs.push_back(std::make_pair(i, j));
                }
            }
        }
    }
}
