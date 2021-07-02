#ifndef RK4_H
#define RK4_H

namespace QUtil {

    namespace algorithm {

        /// Runge-Kutta 4 order method
        /// \tparam T state must impl CopyTo and Accumulate func
        /// \tparam Td derived state must impl CopyTo func
        /// \param func need to be a function pointer rather than lambda
        /// \param state start state
        /// \param stateDerived derived state array
        /// \param tmp_s tmp state
        /// \param tmp_td tmp dervied
        /// \param dt delta t
        template<class T, class Td>
        void RK4(void (*func)(T &, Td &), T &state, Td stateDerived[], T &tmp_s, Td &tmp_td, double dt) {
            double t_arr[4]{dt, dt / 2, dt / 2, dt};
            double w_arr[4]{dt / 6, dt / 3, dt / 3, dt / 6};
            state.CopyTo(tmp_s);
            for (int i = 0; i < 4; ++i) {
                func(tmp_s, stateDerived[i]);
                if (i == 4)
                    break;
                stateDerived[i].CopyTo(tmp_td);
                state.CopyTo(tmp_s);
                tmp_s.Accumulate(tmp_td, t_arr[i]);
            }
            for (int i = 0; i < 4; ++i) {
                state.Accumulate(stateDerived[i], w_arr[i]);
            }
        }
    }

}
#endif