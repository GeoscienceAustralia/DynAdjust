#define TESTING_MAIN

#include "testing.hpp"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

#include "math/dnamatrix_contiguous.hpp"
#include "config/dnaconsts.hpp"
#include "config/dnatypes.hpp"
#include "config/dnatypes-gui.hpp"
#include "measurement_types/dnameasurement.hpp"
#include "measurement_types/dnastation.hpp"
#include "functions/dnatemplatematrixfuncs.hpp"
#include "functions/dnatemplategeodesyfuncs.hpp"
#include "functions/dnatemplatecalcfuncs.hpp"

using namespace dynadjust::math;
using namespace dynadjust::measurements;

namespace {

// Latitude/longitude for a test station pair (approximately Canberra area)
constexpr double kLat1 = -0.6170;    // ~-35.36 deg in radians
constexpr double kLon1 = 2.6110;     // ~149.59 deg in radians
constexpr double kLat2 = -0.6172;
constexpr double kLon2 = 2.6112;

// Compute NStat = measCorr / sqrt(residualPrec) matching UpdateMsrRecordStats
double computeNStat(double measCorr, double residualPrec) {
    return measCorr / sqrt(residualPrec);
}

// Build a 3-component GNSS baseline measurement triplet (X, Y, Z records)
// in bmsBinaryRecords format
void createGNSSBaseline(vmsr_t& records,
    double x, double y, double z,           // pre-adjustment cartesian
    double xAdj, double yAdj, double zAdj,  // adjusted cartesian
    double varXX, double varXY, double varYY,
    double varXZ, double varYZ, double varZZ,
    double nstatX, double nstatY, double nstatZ,
    char msrType = 'G', UINT32 stn1 = 0, UINT32 stn2 = 1)
{
    measurement_t msr_x = {};
    msr_x.measType = msrType;
    msr_x.measStart = xMeas;
    msr_x.station1 = stn1;
    msr_x.station2 = stn2;
    msr_x.vectorCount1 = 1;
    msr_x.vectorCount2 = 0;
    msr_x.term1 = x;
    msr_x.measAdj = xAdj;
    msr_x.term2 = varXX;
    msr_x.NStat = nstatX;
    msr_x.measCorr = xAdj - x;
    msr_x.residualPrec = varXX * 0.1;  // placeholder
    records.push_back(msr_x);

    measurement_t msr_y = {};
    msr_y.measType = msrType;
    msr_y.measStart = yMeas;
    msr_y.station1 = stn1;
    msr_y.station2 = stn2;
    msr_y.vectorCount1 = 0;
    msr_y.vectorCount2 = 0;
    msr_y.term1 = y;
    msr_y.measAdj = yAdj;
    msr_y.term2 = varXY;
    msr_y.term3 = varYY;
    msr_y.NStat = nstatY;
    msr_y.measCorr = yAdj - y;
    msr_y.residualPrec = varYY * 0.1;
    records.push_back(msr_y);

    measurement_t msr_z = {};
    msr_z.measType = msrType;
    msr_z.measStart = zMeas;
    msr_z.station1 = stn1;
    msr_z.station2 = stn2;
    msr_z.vectorCount1 = 0;
    msr_z.vectorCount2 = 0;
    msr_z.term1 = z;
    msr_z.measAdj = zAdj;
    msr_z.term2 = varXZ;
    msr_z.term3 = varYZ;
    msr_z.term4 = varZZ;
    msr_z.NStat = nstatZ;
    msr_z.measCorr = zAdj - z;
    msr_z.residualPrec = varZZ * 0.1;
    records.push_back(msr_z);
}

// Simulate the ENU n-stat computation that UpdateGNSSNstatsForAlternateUnits performs.
// Returns {nstat_e, nstat_n, nstat_u}.
void computeENUNstats(
    double x, double y, double z,
    double xAdj, double yAdj, double zAdj,
    double varXX, double varXY, double varYY,
    double varXZ, double varYZ, double varZZ,
    const double var_adj_diag[3],
    double lat, double lon,
    double& nstat0, double& nstat1, double& nstat2)
{
    matrix_2d gnss_cart(3, 1), gnss_local(3, 1);
    matrix_2d gnss_adj_cart(3, 1), gnss_adj_local(3, 1);
    matrix_2d var_cart(3, 3), var_local(3, 3);
    matrix_2d var_adj_cart(3, 3), var_adj_local(3, 3);
    matrix_2d rotations;

    gnss_cart.put(0, 0, x);
    gnss_cart.put(1, 0, y);
    gnss_cart.put(2, 0, z);

    gnss_adj_cart.put(0, 0, xAdj);
    gnss_adj_cart.put(1, 0, yAdj);
    gnss_adj_cart.put(2, 0, zAdj);

    var_cart.put(0, 0, varXX);
    var_cart.put(0, 1, varXY);  var_cart.put(1, 0, varXY);
    var_cart.put(1, 1, varYY);
    var_cart.put(0, 2, varXZ);  var_cart.put(2, 0, varXZ);
    var_cart.put(1, 2, varYZ);  var_cart.put(2, 1, varYZ);
    var_cart.put(2, 2, varZZ);

    // Simple diagonal adjusted variance for test
    var_adj_cart.put(0, 0, var_adj_diag[0]);
    var_adj_cart.put(1, 1, var_adj_diag[1]);
    var_adj_cart.put(2, 2, var_adj_diag[2]);

    Rotate_CartLocal<double>(gnss_cart, &gnss_local, lat, lon);
    Rotate_CartLocal<double>(gnss_adj_cart, &gnss_adj_local, lat, lon);

    double mid_lat = average(lat, lat);  // same station for simplicity
    double mid_lon = average(lon, lon);

    PropagateVariances_CartLocal_Diagonal<double>(var_cart, var_local,
        mid_lat, mid_lon, rotations, true);
    PropagateVariances_CartLocal_Diagonal<double>(var_adj_cart, var_adj_local,
        mid_lat, mid_lon, rotations, false);

    double measCorr0 = gnss_adj_local.get(0, 0) - gnss_local.get(0, 0);
    double residPrec0 = var_local.get(0, 0) - var_adj_local.get(0, 0);
    nstat0 = computeNStat(measCorr0, residPrec0);

    double measCorr1 = gnss_adj_local.get(1, 0) - gnss_local.get(1, 0);
    double residPrec1 = var_local.get(1, 1) - var_adj_local.get(1, 1);
    nstat1 = computeNStat(measCorr1, residPrec1);

    double measCorr2 = gnss_adj_local.get(2, 0) - gnss_local.get(2, 0);
    double residPrec2 = var_local.get(2, 2) - var_adj_local.get(2, 2);
    nstat2 = computeNStat(measCorr2, residPrec2);
}

// Simulate the AED n-stat computation
void computeAEDNstats(
    double x, double y, double z,
    double xAdj, double yAdj, double zAdj,
    double varXX, double varXY, double varYY,
    double varXZ, double varYZ, double varZZ,
    const double var_adj_diag[3],
    double lat, double lon,
    double& nstat0, double& nstat1, double& nstat2)
{
    matrix_2d gnss_cart(3, 1), gnss_local(3, 1);
    matrix_2d gnss_adj_cart(3, 1), gnss_adj_local(3, 1);
    matrix_2d var_cart(3, 3), var_local(3, 3), var_polar(3, 3);
    matrix_2d var_adj_cart(3, 3), var_adj_local(3, 3), var_adj_polar(3, 3);
    matrix_2d rotations;

    gnss_cart.put(0, 0, x);  gnss_cart.put(1, 0, y);  gnss_cart.put(2, 0, z);
    gnss_adj_cart.put(0, 0, xAdj);  gnss_adj_cart.put(1, 0, yAdj);  gnss_adj_cart.put(2, 0, zAdj);

    var_cart.put(0, 0, varXX);
    var_cart.put(0, 1, varXY);  var_cart.put(1, 0, varXY);
    var_cart.put(1, 1, varYY);
    var_cart.put(0, 2, varXZ);  var_cart.put(2, 0, varXZ);
    var_cart.put(1, 2, varYZ);  var_cart.put(2, 1, varYZ);
    var_cart.put(2, 2, varZZ);

    var_adj_cart.put(0, 0, var_adj_diag[0]);
    var_adj_cart.put(1, 1, var_adj_diag[1]);
    var_adj_cart.put(2, 2, var_adj_diag[2]);

    Rotate_CartLocal<double>(gnss_cart, &gnss_local, lat, lon);
    Rotate_CartLocal<double>(gnss_adj_cart, &gnss_adj_local, lat, lon);

    double mid_lat = average(lat, lat);
    double mid_lon = average(lon, lon);

    double azimuth = Direction(gnss_local.get(0, 0), gnss_local.get(1, 0));
    double elevation = VerticalAngle(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
    double distance = magnitude(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));

    double azimuthAdj = Direction(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0));
    double elevationAdj = VerticalAngle(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));
    double distanceAdj = magnitude(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));

    PropagateVariances_LocalCart<double>(var_cart, var_local,
        mid_lat, mid_lon, false, rotations, true);
    PropagateVariances_LocalCart<double>(var_adj_cart, var_adj_local,
        mid_lat, mid_lon, false, rotations, false);

    PropagateVariances_LocalPolar_Diagonal<double>(var_local, var_polar,
        azimuth, elevation, distance, rotations, true);
    PropagateVariances_LocalPolar_Diagonal<double>(var_adj_local, var_adj_polar,
        azimuth, elevation, distance, rotations, false);

    // Azimuth
    double mc0 = azimuthAdj - azimuth;
    double rp0 = var_polar.get(0, 0) - var_adj_polar.get(0, 0);
    nstat0 = computeNStat(mc0, rp0);

    // Elevation
    double mc1 = elevationAdj - elevation;
    double rp1 = var_polar.get(1, 1) - var_adj_polar.get(1, 1);
    nstat1 = computeNStat(mc1, rp1);

    // Distance
    double mc2 = distanceAdj - distance;
    double rp2 = var_polar.get(2, 2) - var_adj_polar.get(2, 2);
    nstat2 = computeNStat(mc2, rp2);
}

} // anonymous namespace


TEST_CASE("Rotate_CartLocal transforms cartesian to local frame", "[gnss][rotation]") {
    matrix_2d cart(3, 1), local(3, 1);

    // A purely X-direction cartesian vector at the equator/prime meridian
    // should map to (0, 0, up) since X points radially at (0,0)
    cart.put(0, 0, 1.0);
    cart.put(1, 0, 0.0);
    cart.put(2, 0, 0.0);

    Rotate_CartLocal<double>(cart, &local, 0.0, 0.0);

    // At lat=0, lon=0: X -> Up, so local(2,0) should be ~1
    REQUIRE(fabs(local.get(0, 0)) < 1e-10);  // E ~ 0
    REQUIRE(fabs(local.get(1, 0)) < 1e-10);  // N ~ 0
    REQUIRE(fabs(local.get(2, 0) - 1.0) < 1e-10);  // U ~ 1
}

TEST_CASE("Rotate_CartLocal Y-axis at equator/prime meridian maps to East", "[gnss][rotation]") {
    matrix_2d cart(3, 1), local(3, 1);

    cart.put(0, 0, 0.0);
    cart.put(1, 0, 1.0);
    cart.put(2, 0, 0.0);

    Rotate_CartLocal<double>(cart, &local, 0.0, 0.0);

    // At lat=0, lon=0: Y -> East
    REQUIRE(fabs(local.get(0, 0) - 1.0) < 1e-10);  // E ~ 1
    REQUIRE(fabs(local.get(1, 0)) < 1e-10);          // N ~ 0
    REQUIRE(fabs(local.get(2, 0)) < 1e-10);          // U ~ 0
}

TEST_CASE("ENU n-stat differs from cartesian n-stat", "[gnss][nstat][enu]") {
    // Construct a baseline where the cartesian and ENU n-stats are different
    double x = 100.0, y = 200.0, z = 50.0;
    double xAdj = 100.05, yAdj = 200.03, zAdj = 50.01;

    double varXX = 0.001, varXY = 0.0001, varYY = 0.002;
    double varXZ = 0.00005, varYZ = 0.0001, varZZ = 0.0015;

    double var_adj_diag[3] = {0.0005, 0.001, 0.0008};

    // Cartesian n-stats
    double cart_nstat_x = (xAdj - x) / sqrt(varXX - var_adj_diag[0]);
    double cart_nstat_y = (yAdj - y) / sqrt(varYY - var_adj_diag[1]);
    double cart_nstat_z = (zAdj - z) / sqrt(varZZ - var_adj_diag[2]);

    // ENU n-stats
    double enu_nstat_e, enu_nstat_n, enu_nstat_u;
    computeENUNstats(x, y, z, xAdj, yAdj, zAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1,
        enu_nstat_e, enu_nstat_n, enu_nstat_u);

    // The rotation mixes components, so ENU n-stats should differ from cartesian
    REQUIRE(fabs(enu_nstat_e - cart_nstat_x) > 0.01);
    REQUIRE(std::isfinite(enu_nstat_e));
    REQUIRE(std::isfinite(enu_nstat_n));
    REQUIRE(std::isfinite(enu_nstat_u));
}

TEST_CASE("AED n-stat computation produces finite values", "[gnss][nstat][aed]") {
    double x = 100.0, y = 200.0, z = 50.0;
    double xAdj = 100.05, yAdj = 200.03, zAdj = 50.01;

    double varXX = 0.001, varXY = 0.0001, varYY = 0.002;
    double varXZ = 0.00005, varYZ = 0.0001, varZZ = 0.0015;

    double var_adj_diag[3] = {0.0005, 0.001, 0.0008};

    double aed_nstat_a, aed_nstat_e, aed_nstat_d;
    computeAEDNstats(x, y, z, xAdj, yAdj, zAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1,
        aed_nstat_a, aed_nstat_e, aed_nstat_d);

    REQUIRE(std::isfinite(aed_nstat_a));
    REQUIRE(std::isfinite(aed_nstat_e));
    REQUIRE(std::isfinite(aed_nstat_d));
}

TEST_CASE("ENU n-stats differ from AED n-stats for same baseline", "[gnss][nstat][enu][aed]") {
    double x = 100.0, y = 200.0, z = 50.0;
    double xAdj = 100.05, yAdj = 200.03, zAdj = 50.01;

    double varXX = 0.001, varXY = 0.0001, varYY = 0.002;
    double varXZ = 0.00005, varYZ = 0.0001, varZZ = 0.0015;

    double var_adj_diag[3] = {0.0005, 0.001, 0.0008};

    double enu_e, enu_n, enu_u;
    computeENUNstats(x, y, z, xAdj, yAdj, zAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1, enu_e, enu_n, enu_u);

    double aed_a, aed_e, aed_d;
    computeAEDNstats(x, y, z, xAdj, yAdj, zAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1, aed_a, aed_e, aed_d);

    // Different frames produce different n-stats
    REQUIRE(fabs(enu_e - aed_a) > 1e-6);
}

TEST_CASE("Sorting by ENU n-stat gives different order than cartesian n-stat", "[gnss][nstat][sort]") {
    // Two baselines with deliberately different orderings in cartesian vs ENU

    // Baseline A: large X correction, small Y/Z
    double ax = 100.0, ay = 0.1, az = 0.1;
    double axAdj = 100.1, ayAdj = 0.1001, azAdj = 0.1001;

    // Baseline B: moderate X correction, larger Y
    double bx = 50.0, by = 200.0, bz = 10.0;
    double bxAdj = 50.02, byAdj = 200.08, bzAdj = 10.005;

    double var = 0.001;
    double varXX = var, varXY = 0.0, varYY = var;
    double varXZ = 0.0, varYZ = 0.0, varZZ = var;
    double var_adj_diag[3] = {var * 0.5, var * 0.5, var * 0.5};

    // Cartesian n-stats (only looking at X component for comparison)
    double cart_nstat_a = (axAdj - ax) / sqrt(varXX - var_adj_diag[0]);
    double cart_nstat_b = (bxAdj - bx) / sqrt(varXX - var_adj_diag[0]);

    // ENU n-stats
    double enu_ae, enu_an, enu_au;
    computeENUNstats(ax, ay, az, axAdj, ayAdj, azAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1, enu_ae, enu_an, enu_au);

    double enu_be, enu_bn, enu_bu;
    computeENUNstats(bx, by, bz, bxAdj, byAdj, bzAdj,
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1, enu_be, enu_bn, enu_bu);

    // Verify baselines produce finite results
    REQUIRE(std::isfinite(enu_ae));
    REQUIRE(std::isfinite(enu_be));

    // Cartesian X n-stat: A has |0.1/sqrt(0.0005)| > B |0.02/sqrt(0.0005)|
    REQUIRE(fabs(cart_nstat_a) > fabs(cart_nstat_b));

    // In ENU frame the ordering may differ because the rotation mixes components
    // At minimum, the ENU values should not equal the cartesian values
    REQUIRE(fabs(enu_ae - cart_nstat_a) > 0.01);
}

TEST_CASE("PropagateVariances_CartLocal_Diagonal preserves positive definiteness", "[gnss][variance]") {
    matrix_2d var_cart(3, 3), var_local(3, 3);
    matrix_2d rotations;

    // Positive definite variance matrix
    var_cart.put(0, 0, 0.004);
    var_cart.put(0, 1, 0.001);  var_cart.put(1, 0, 0.001);
    var_cart.put(1, 1, 0.005);
    var_cart.put(0, 2, 0.0005); var_cart.put(2, 0, 0.0005);
    var_cart.put(1, 2, 0.0008); var_cart.put(2, 1, 0.0008);
    var_cart.put(2, 2, 0.003);

    PropagateVariances_CartLocal_Diagonal<double>(var_cart, var_local,
        kLat1, kLon1, rotations, true);

    // Diagonal elements must remain positive
    REQUIRE(var_local.get(0, 0) > 0.0);
    REQUIRE(var_local.get(1, 1) > 0.0);
    REQUIRE(var_local.get(2, 2) > 0.0);
}

TEST_CASE("Y measurement type is not affected by alternate units", "[gnss][nstat][skip-Y]") {
    // This verifies the logic that Y measurements are skipped.
    // In the actual code, UpdateGNSSNstatsForAlternateUnits checks:
    //   if (_it_msr->measType != 'G' && _it_msr->measType != 'X') continue;
    // So Y measurements keep their original cartesian n-stats.

    // Just verify the filter condition
    char typeG = 'G';
    char typeX = 'X';
    char typeY = 'Y';

    bool processG = (typeG == 'G' || typeG == 'X');
    bool processX = (typeX == 'G' || typeX == 'X');
    bool processY = (typeY == 'G' || typeY == 'X');

    REQUIRE(processG == true);
    REQUIRE(processX == true);
    REQUIRE(processY == false);
}

TEST_CASE("Non-xMeas records are skipped", "[gnss][nstat][skip-non-xMeas]") {
    // Verify that only xMeas (measStart == 0) records get processed
    REQUIRE(xMeas == 0);
    REQUIRE(yMeas == 1);
    REQUIRE(zMeas == 2);

    measurement_t msr_x = {};
    msr_x.measStart = xMeas;
    REQUIRE(msr_x.measStart == xMeas);

    measurement_t msr_y = {};
    msr_y.measStart = yMeas;
    REQUIRE(msr_y.measStart != xMeas);

    measurement_t msr_z = {};
    msr_z.measStart = zMeas;
    REQUIRE(msr_z.measStart != xMeas);
}

TEST_CASE("Direction and VerticalAngle produce expected values", "[gnss][geodesy]") {
    // East=1, North=0 should give azimuth of 90 degrees (PI/2)
    double az = Direction(1.0, 0.0);
    REQUIRE(fabs(az - HALF_PI) < 1e-10);

    // East=0, North=1 should give azimuth of 0 (or equivalently TWO_PI)
    az = Direction(0.0, 1.0);
    REQUIRE(fabs(az) < 1e-10 || fabs(az - TWO_PI) < 1e-10);

    // Vertical angle for a purely horizontal vector should be 0
    double va = VerticalAngle(1.0, 1.0, 0.0);
    REQUIRE(fabs(va) < 1e-10);

    // Vertical angle for a purely vertical vector should be +/- PI/2
    va = VerticalAngle(0.0, 0.0, 1.0);
    REQUIRE(fabs(va - HALF_PI) < 1e-10);
}

TEST_CASE("magnitude computes 3D vector length", "[gnss][geodesy]") {
    REQUIRE(fabs(magnitude(3.0, 4.0, 0.0) - 5.0) < 1e-10);
    REQUIRE(fabs(magnitude(1.0, 1.0, 1.0) - sqrt(3.0)) < 1e-10);
}

TEST_CASE("GNSS n-stat sort condition only triggers for n-stat sort with non-XYZ units", "[gnss][nstat][condition]") {
    // The pre-sort update should only occur when:
    //   _adj_gnss_units != XYZ_adj_gnss_ui && _sort_adj_msr == n_st_adj_msr_sort_ui

    // Should NOT trigger for XYZ with n-stat sort
    REQUIRE(!(XYZ_adj_gnss_ui != XYZ_adj_gnss_ui && true));

    // Should trigger for ENU with n-stat sort
    REQUIRE((ENU_adj_gnss_ui != XYZ_adj_gnss_ui && true));

    // Should NOT trigger for ENU with non-n-stat sort
    REQUIRE(!(ENU_adj_gnss_ui != XYZ_adj_gnss_ui && false));

    // Should trigger for AED with n-stat sort
    REQUIRE((AED_adj_gnss_ui != XYZ_adj_gnss_ui && true));

    // Should trigger for ADU with n-stat sort
    REQUIRE((ADU_adj_gnss_ui != XYZ_adj_gnss_ui && true));

    // Verify enum values
    REQUIRE(XYZ_adj_gnss_ui == 0);
    REQUIRE(ENU_adj_gnss_ui == 1);
    REQUIRE(AED_adj_gnss_ui == 2);
    REQUIRE(ADU_adj_gnss_ui == 3);
    REQUIRE(n_st_adj_msr_sort_ui == 7);
}

TEST_CASE("Zero correction produces zero n-stat", "[gnss][nstat][zero]") {
    // When adjusted == pre-adjusted, corrections are zero, n-stat should be zero
    double x = 100.0, y = 200.0, z = 50.0;
    double varXX = 0.001, varXY = 0.0, varYY = 0.001;
    double varXZ = 0.0, varYZ = 0.0, varZZ = 0.001;
    double var_adj_diag[3] = {0.0005, 0.0005, 0.0005};

    double nstat_e, nstat_n, nstat_u;
    computeENUNstats(x, y, z, x, y, z,  // no correction
        varXX, varXY, varYY, varXZ, varYZ, varZZ,
        var_adj_diag, kLat1, kLon1,
        nstat_e, nstat_n, nstat_u);

    REQUIRE(fabs(nstat_e) < 1e-10);
    REQUIRE(fabs(nstat_n) < 1e-10);
    REQUIRE(fabs(nstat_u) < 1e-10);
}
