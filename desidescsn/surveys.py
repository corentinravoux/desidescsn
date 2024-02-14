# DESI
# Hypothesis:
# - BGS only
# - No very bright limit cut
# - No rfiber cut
# - No differenciation between BGS Bright and Faint
# - Only one Fiber assignement efficiency both for Bright and Faint
# - Minimal fiber efficiency (conservative choice)
# - No color cut (including stellar rejection)

magnitude_cut_DESIBGS = 20.175
band_DESIBGS = "r"
fiber_assignement_efficiency_DESIBGS = 0.8
redshift_success_rate_DESIBGS = 0.95

# DESI2
# Hypothesis:
# - BGS only
# - zfiber = z (very strong and can be removed)
# - No color cut (including stellar rejection)

magnitude_cut_DESI2 = 21.6
band_DESI2 = "zfiber"
fiber_assignement_efficiency_DESI2 = 0.8
redshift_success_rate_DESI2 = 0.95


# 4MOST CRS
# Hypothesis:
# - BG only
# - J_VHS = J_ROMAN
# - No W1 color cut
# - No fiber assignement efficiency loss

magnitude_cut_CRS_1 = 16
magnitude_cut_CRS_2 = 18
color_cut_CRS = 0.1
band_CRS_1 = "J"
band_CRS_2 = "K"

fiber_assignement_efficiency_CRS = 1.0
redshift_success_rate_CRS = 0.95


# 4MOST 4HS
# Hypothesis:
# - J_VHS = J_ROMAN
# - K_VHS = K_ROMAN
# - Only one color cut
# - No PV sub-sample included
# - No fiber assignement efficiency loss


magnitude_cut_4HS = 18
color_cut_4HS = 0.45
band_4HS_1 = "J"
band_4HS_2 = "K"

fiber_assignement_efficiency_4HS = 1.0
redshift_success_rate_4HS = 0.98


def mask_magnitude_DESIBGS(
    file_sn,
    hashing_table,
):
    key_r = hashing_table[band_DESIBGS]
    mask_magnitude = file_sn[key_r] < magnitude_cut_DESIBGS

    return mask_magnitude


def host_efficiency_DESIBGS():
    return fiber_assignement_efficiency_DESIBGS * redshift_success_rate_DESIBGS


def mask_magnitude_DESI2(
    file_sn,
    hashing_table,
):
    key_band = hashing_table[band_DESI2]
    mask_magnitude = file_sn[key_band] < magnitude_cut_DESI2

    return mask_magnitude


def host_efficiency_DESI2():
    return fiber_assignement_efficiency_DESI2 * redshift_success_rate_DESI2


def mask_magnitude_CRS(
    file_sn,
    hashing_table,
):
    key_band_1 = hashing_table[band_CRS_1]
    key_band_2 = hashing_table[band_CRS_2]

    color_to_cut = file_sn[key_band_1] - file_sn[key_band_2]

    mask_magnitude = file_sn[key_band_1] > magnitude_cut_CRS_1
    mask_magnitude &= file_sn[key_band_1] < magnitude_cut_CRS_2
    mask_magnitude &= color_to_cut > color_cut_CRS

    return mask_magnitude


def host_efficiency_CRS():
    return fiber_assignement_efficiency_CRS * redshift_success_rate_CRS


def mask_magnitude_4HS(
    file_sn,
    hashing_table,
):
    key_band_1 = hashing_table[band_4HS_1]
    key_band_2 = hashing_table[band_4HS_2]

    color_to_cut = file_sn[key_band_1] - file_sn[key_band_2]
    mask_magnitude = file_sn[key_band_1] < magnitude_cut_4HS
    mask_magnitude &= color_to_cut < color_cut_4HS

    return mask_magnitude


def host_efficiency_4HS():
    return fiber_assignement_efficiency_4HS * redshift_success_rate_4HS
