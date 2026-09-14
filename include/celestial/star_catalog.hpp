// star_catalog.hpp — Yale Bright Star Catalogue (BSC5).
//
// 5080 stars to V = 6.0, generated into star_catalog_data.cpp by
// tools/analysis/make_catalog.py. Identifiers are HR (Harvard Revised) numbers.
//
// The 53.5 x 35.1 deg field of the paper's camera covers ~4.5% of the sky, so
// the full V=6.0 list puts ~230 stars in frame. That is far more than any real
// detector will pull out of a 100 ms exposure, so set the DETECTION limit via
// ErrorModel::mag_limit rather than by truncating the catalogue: it keeps the
// astrometry honest and lets you sweep detector sensitivity as a parameter.
//
// Rough guide for the paper's f/1.4 6 mm optics at 10 Hz: V ~ 4.0 gives about
// 25 stars in frame, matching the 27 tracked in their Figure 2.

#pragma once

#include <vector>

#include "celestial/types.hpp"

namespace celestial {

namespace detail {
/// Raw generated row. See tools/analysis/make_catalog.py.
struct CatalogRaw {
  int hr;
  double ra_deg;
  double dec_deg;
  double pm_ra_arcsec;   ///< dRA/dt, arcsec/yr. NOT mu_alpha*.
  double pm_dec_arcsec;  ///< dDec/dt, arcsec/yr
  double vmag;
};
extern const CatalogRaw kBsc5[];
extern const int kBsc5Count;
}  // namespace detail

/// The full built-in catalogue, sorted by magnitude (brightest first).
const std::vector<CatalogStar>& brightStarCatalog();

/// Subset brighter than a given visual magnitude. Because the catalogue is
/// magnitude-sorted this is a prefix, so it is cheap.
std::vector<CatalogStar> catalogBrighterThan(double vmag_limit);

/// O(1) lookup by catalogue id. Returns nullptr if unknown.
const CatalogStar* findStar(int id);

/// Common name for a catalogue id, or nullptr. Only the ~65 navigational
/// stars are named; callers should fall back to "HR nnnn".
const char* starName(int id);

}  // namespace celestial
