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
        
        m = static_cast<std::int32_t>(periodiclen / SL) - 1;
        mesh_size = static_cast<double>(periodiclen) / m;
        
        BOOST_ASSERT(m > 2);
        BOOST_ASSERT(mesh_size > SL);
        
        number_of_mesh = m * m * m;
        count.resize(number_of_mesh);
        indexes.resize(number_of_mesh);
    }

    void MeshList::make_pair(SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs)
    {
        pairs.clear();
        
        auto const pn = atoms.size();

        std::vector<std::int32_t> particle_position(pn);
        std::vector<std::int32_t> pointer(number_of_mesh);
        std::fill(particle_position.begin(), particle_position.end(), 0);
        std::fill(count.begin(), count.end(), 0);
        std::fill(pointer.begin(), pointer.end(), 0);

        auto const im = 1.0 / mesh_size;
        for (auto i = 0; i < pn; i++) {
            auto ix = static_cast<std::int32_t>(atoms[i].r[0] * im);
            auto iy = static_cast<std::int32_t>(atoms[i].r[1] * im);
            auto iz = static_cast<std::int32_t>(atoms[i].r[2] * im);
            
            if (ix < 0) {
                ix += m;
            }
            else if (ix >= m) {
                ix -= m;
            }

            if (iy < 0) {
                iy += m;
            }
            else if (iy >= m) {
                iy -= m;
            }
            if (iz < 0) {
                iz += m;
            }
            else if (iz >= m) {
                iz -= m;
            }

            auto const index = ix + iy * m + iz * m * m;
            
            BOOST_ASSERT(index >= 0);
            BOOST_ASSERT(index < number_of_mesh);
            
            count[index]++;
            particle_position[i] = index;
        }
        
        indexes[0] = 0;
        auto sum = 0;
        
        for (auto i = 0; i < number_of_mesh - 1; i++) {
            sum += count[i];
            indexes[i + 1] = sum;
        }
        
        for (auto i = 0; i < pn; i++) {
            auto const pos = particle_position[i];
            auto const j = indexes[pos] + pointer[pos];
            sorted_buffer[j] = i;
            ++pointer[pos];
        }
        
        for (auto i = 0; i < number_of_mesh; i++) {
            search(i, atoms, pairs);
        }
    }

    void MeshList::search_other(std::int32_t id, std::int32_t ix, std::int32_t iy, std::int32_t iz, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs)
    {
        if (ix < 0) {
            ix += m;
        }
        else if (ix >= m) {
            ix -= m;
        }
        
        if (iy < 0) {
            iy += m;
        }
        else if (iy >= m) {
            iy -= m;
        }
        
        if (iz < 0) {
            iz += m;
        }
        else if (iz >= m) {
            iz -= m;
        }
        
        auto const id2 = ix + iy * m + iz * m * m;

        for (auto k = indexes[id]; k < indexes[id] + count[id]; k++) {
            for (auto m = indexes[id2]; m < indexes[id2] + count[id2]; m++) {
                auto const i = sorted_buffer[k];
                auto const j = sorted_buffer[m];

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
        auto const ix = id % m;
        auto const iy = (id / m) % m;
        auto const iz = (id / m / m);

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
        auto const si = indexes[id];
        auto const n = count[id];
        for (auto k = si; k < si + n - 1; k++) {
            for (auto m = k + 1; m < si + n; m++) {
                auto const i = sorted_buffer[k];
                auto const j = sorted_buffer[m];

                Eigen::Vector4d d = atoms[j].r - atoms[i].r;

                SystemParam::adjust_periodic(d, periodiclen_);

                if (d.squaredNorm() <= SystemParam::ML2) {
                    pairs.push_back(std::make_pair(i, j));
                }
            }
        }
    }
}
