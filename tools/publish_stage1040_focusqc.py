from pathlib import Path
import pandas as pd
P=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review");F=P/"single-cell.html";s=F.read_text(encoding="utf8")
img="assets/stage1038_threegene_focusQC/STAGE1038_THREE_GENE_NATIVE_XY_Z_SHARPNESS.png"
assert (P/img).is_file()
A=pd.read_csv(P/"data"/"stage1038_n5_native_spot_sharpness_fwhm_summary.csv")
by=A.set_index("gene")
assert int(by.loc["VGAT","n"])>=80
t=f"""<section id="stage1038-vgat-focus"><h2>Does VGAT really appear more defocused? Native-voxel 3D spot-width check</h2>
<p class="lead"><b>At high native resolution the detected VGAT puncta are not systematically wider in this sample.</b> We randomly sampled detected source spots and retained ones with local peak signal-to-background &gt;4 to compare true bright-spot widths in the original 0.23µm XY and 0.42µm Z voxels. Measured median lateral punctum FWHM: VGAT <b>{by.loc['VGAT','xy_FWHM_median_um']:.2f}µm</b> (n={int(by.loc['VGAT','n'])}), VGLUT2 {by.loc['VGLUT2','xy_FWHM_median_um']:.2f}µm (n={int(by.loc['VGLUT2','n'])}), Snap25 {by.loc['Snap25','xy_FWHM_median_um']:.2f}µm (n={int(by.loc['Snap25','n'])}). Axial median widths: {by.loc['VGAT','z_FWHM_median_um']:.2f}, {by.loc['VGLUT2','z_FWHM_median_um']:.2f}, and {by.loc['Snap25','z_FWHM_median_um']:.2f}µm respectively. This argues against a universal VGAT-specific loss of focus <em>among strong detected puncta in one S500 source image</em>. It does not test dim uncalled VGAT puncta or all Z planes; those may remain blurred and underdetected.</p>
<div class="cards"><article class="card"><h3>VGAT, VGLUT2, Snap25 original N5 spot sharpness, XY and Z</h3><a href="{img}"><img loading="lazy" src="{img}"></a><p>Random detector-listed native spots with local high-SNR, widths estimated by contiguous half-height local raw image intensity (not an optical PSF fit). <a href="data/stage1038_n5_native_spot_sharpness_fwhm_summary.csv">Summary values and sample n</a> · <a href="data/stage1038_native_focusQC_authority.json">Selection limits</a>.</p></article></div></section>"""
if "stage1038-vgat-focus" not in s:
 assert '<section id="stage1036-full-raw-spots">' in s
 s=s.replace('<section id="stage1036-full-raw-spots">',t+'<section id="stage1036-full-raw-spots">',1)
F.write_text(s,encoding="utf8")
print("STAGE1040_FOCUS_PAGE_READY",len(s),flush=True)
