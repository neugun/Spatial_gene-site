from pathlib import Path
P=Path(r'G:\Spatial_gene_site_publish\perilc-map6-review')
F=P/'single-cell.html';h=F.read_text(encoding='utf8')
def c(title,p,desc):
 assert (P/p).exists(),p
 return '<article class="card"><h3>'+title+'</h3><a href="'+p+'"><img loading="lazy" src="'+p+'"></a><p>'+desc+'</p></article>'
a='assets/stage1043_strict_transmitter_gate/'
z='assets/stage1044_axial_spot_rescue/'
cards=[
 c('Entire 3-section anatomy after exact gate',a+'STAGE1043_USER_GATE_ALL_SECTIONS_SPATIAL.png','True Stage631 entire physical XY. Grey cells excluded, VGAT-only blue, VGLUT2-only orange, double-positive purple. S500/S530/S560 included.'),
 c('Exact-gate cells on ORIGINAL all-cell t-SNE',a+'STAGE1043_USER_GATE_ORIGINAL_STAGE995_TSNE.png','Same 71,950-cell Stage995 label-free embedding with user filter overlaid after fitting; secondary panel adds existing Allen reference-neuron status. Not retrained here.'),
 c('Threshold robustness: dual fractions and sample size',a+'STAGE1045_GATE_THRESHOLD_SENSITIVITY.png','Snap25>5 and frozen NE/ChAT excluded for all 42 thresholds. Moving cutoffs changes observed dual fractions, but does not validate subtypes.'),
 c('Do VGAT unmatched candidates appear in upper/lower planes?',z+'STAGE1044_UNMATCHED_VGAT_UP_DOWN_Z_RESCUE.png','Region 1: 86 native raw VGAT candidate peaks. Full original RS-FISH spot list searched at same XY (within 3 original pixels) through ±4/8/16/32/64 z voxels. 42 match original 3D; of 44 otherwise unmatched, zero have a spot in adjacent ±32 layers; one very distant XY-similar detection only within ±64. This does not constitute the same molecule.')
]
body='<section id="stage1043-user-strict-gate"><h2>Strict neurotransmitter-neuron-only sensitivity and VGAT up/down-Z test</h2>'
body+='<p class="lead"><b>Actual gate on all 71,950 cells:</b> corrected Snap25 &gt;5 AND frozen broad group neither NE nor ChAT AND (VGAT/Slc32a1 &gt;10 OR VGLUT2/Slc17a6 &gt;5). This retains <b>31,115 cells</b> (43.25%), split into VGAT-only <b>12,886</b>, VGLUT2-only <b>7,257</b>, and <b>10,972 (35.26%) dual-marker</b>. The gate makes a useful operational transmitter-positive population, <b>but cannot make double-positives disappear</b>; selecting only marker-positive cells increases their fraction among remaining cells.</p>'
body+='<p class="lead"><b>Original taxonomy caveat:</b> 24,759 retained are old inhibitory-broad cells and 6,346 are old excitatory-broad; the remainder are 10 old Other. Among old inhibitory cells 9,890/24,759 (39.95%) are double-positive versus 1,081/6,346 (17.03%) among excitatory. The broad classes themselves were inferred from these genes, not independent ground truth. An extra old reference-neuron intersection retains 30,830 cells with 10,863 duals: the decision barely changes. Excluding NE and ChAT is based on frozen Broad identity, NOT forcing zero Slc6a2/Slc5a7; <a href="data/stage1043_special_marker_exclusion_sensitivity.csv">marker-wise special-class exclusion differs sharply</a>.</p>'
body+='<p class="lead"><b>Z-plane answer:</b> using the complete VGAT source-round 3D spot list in the same source XY field, the 44 unmatched strong image maxima were rechecked across ±4, ±8, ±16, ±32 and ±64 original Z layers (each 0.42µm). <b>None are rescued within ±32 Z planes</b>; one candidate has a distinct XY-near spot at very distant Z within ±64. The broad candidate test expanded the z search and finds 42 original 3D matches instead of the prior 39, because the edge of the original 15-plane interval was excluded previously. These image maxima are still not validated uncalled RNA.</p>'
body+='<p class="lead">If preserving only directly exclusive marker populations (VGAT-only + VGLUT2-only), the sample drops to <b>20,143 cells</b> and discards 35.26% of currently retained marker-positive cells. This may help make a readable two-color visualization but would systematically exclude biologically ambiguous cells, not resolve the problem. Higher VGAT/VGLUT2 cutoffs also trade sample size for lower apparent double rates (e.g. 10/10: 26,685 cells, 19.5% double). No forced dominant-transmitter assignment is made.</p>'
body+='<div class="cards">'+''.join(cards)+'</div>'
body+='<p class="lead"><a href="data/stage1043_exact_gate_counts.csv">Exact numbers for 3 sections and reference sensitivity</a> · <a href="data/stage1043_frozen_broad_versus_actual_marker_call.csv">Old E/I versus actual marker call</a> · <a href="data/stage1045_exact_gate_threshold_sensitivity.csv">All thresholds</a> · <a href="data/stage1045_dual_strength.json">Strong-versus-borderline coexpression</a> · <a href="data/stage1044_near_xy_upperlower_plane_rescue_counts.csv">Up/down Z search results</a>. All come from the same mouse (three sections), not independent animal replication.</p></section>'
if 'stage1043-user-strict-gate' not in h:
 anchor='<section id="stage1041-candidate-recall">'
 assert anchor in h
 h=h.replace(anchor,body+anchor,1)
 h=h.replace('<a href="#stage1036-full-raw-spots">','<a href="#stage1043-user-strict-gate">Strict neuron-only gate / Z rescue</a><a href="#stage1036-full-raw-spots">',1)
F.write_text(h,encoding='utf8')
print('STAGE1047_SITE_READY',len(h),len(cards),flush=True)
