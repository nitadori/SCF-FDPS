// EPJ-Plummer
template <>
void CalcGravityP0(const EPI * ep_i,
                   const PS::S32 n_ip,
                   const EPJ * ep_j,
                   const PS::S32 n_jp,
                   Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    const auto epsinv = 1.0 / EPI::eps;
    nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

// SPJ-P0-Plummer
template <>
void CalcGravityP0(const EPI * ep_i,
                   const PS::S32 n_ip,
                   const SPJ * ep_j,
                   const PS::S32 n_jp,
                   Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    const auto epsinv = 1.0 / EPI::eps;
    nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

// SPJ-P2-Plummer
template <>
void CalcGravityP2(const EPI * ep_i,
                    const PS::S32 n_ip,
                    const SPJ * ep_j,
                    const PS::S32 n_jp,
                    Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    nbody_m256d_quad(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_spj(n_ip, n_jp);
}

// EPJ-Spline
template <>
void CalcGravitySpl(const EPI * ep_i,
                    const PS::S32 n_ip,
                    const EPJ * ep_j,
                    const PS::S32 n_jp,
                    Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    const auto epsinv = 1.0 / EPI::eps;
    nbody_m256d_spline(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

// SPJ-P0-Spline
template <>
void CalcGravityP0Eps0(const EPI * ep_i,
                       const PS::S32 n_ip,
                       const SPJ * ep_j,
                       const PS::S32 n_jp,
                       Force * force){
    const auto eps2 = 0.0;
    const auto epsinv = 0.0;

    nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

// SPJ-P2-Spline
template <>
void CalcGravityP2Eps0(const EPI * ep_i,
                       const PS::S32 n_ip,
                       const SPJ * ep_j,
                       const PS::S32 n_jp,
                       Force * force){
    const auto eps2 = 0.0;
    nbody_m256d_quad(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_spj(n_ip, n_jp);
}
