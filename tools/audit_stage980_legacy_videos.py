"""Map10 Stage980 reconstruct early video standards from verified Stage191/214/215 films."""
from pathlib import Path
import json,subprocess
base=Path(r"G:\PeriLC_current\D_mirror\11_PPTsummary\PeriLC_CHATGPT_RESULTS_20260914")
names=[("S191","191_MAP4_K5_VIDEO","FIG2A_K5_CORONAL_ZOOM_ROTATE.mp4"),("S214","214_MAP4_K5_PYVISTA_MOVIE","FIG2A_K5_SMOOTHBODY_ZOOM_ROTATE.mp4"),("S215","215_MAP4_K5_PYVISTA_MOVIE","K5_SMOOTHBODY_ZOOM_ROTATE.mp4")]
out=[]
for stage,folder,name in names:
 path=base/folder/name
 assert path.is_file()
 p=subprocess.run(["ffprobe","-v","error","-show_entries","format=duration,size:stream=width,height,r_frame_rate,nb_frames,codec_name","-of","json",str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 assert p.returncode==0,(stage,p.stderr)
 obj=json.loads(p.stdout)
 stream=next(s for s in obj["streams"] if s.get("width"))
 shots=sorted(q.name for q in (base/folder).glob("*frame*.png"))
 out.append({"historical_stage":stage,"film":name,"duration_seconds":round(float(obj["format"]["duration"]),3),"file_bytes":int(obj["format"]["size"]),"width":int(stream["width"]),"height":int(stream["height"]),"frame_rate":stream.get("r_frame_rate"),"frames":stream.get("nb_frames"),"available_QC_frame_stills":shots})
report={"stage":980,"historical_video_inventory":out,"current_required_video_gates":["Show full coronal context before any zoom","Smooth camera zoom into selected region; hold view long enough to examine boundaries","True 3D orbital motion for zoomed segmented bodies, not 2D image rotation","Stage631 already x-flipped/y-flipped physical coordinates; no second flip","Use Stage977 Fine26 26-color key consistently in 2D and 3D","Include start, zoom transition, mid-orbit and final QC stills","Verify no clipped tissue, jump cuts, floating artefacts, or disappearing labels"],"status":"HISTORICAL_REFERENCE_VERIFIED_NEW_MOVIE_PENDING"}
dest=Path(r"G:\Spatial_gene_site_publish\perilc-map6-review\data\stage980_legacy_video_reference.json")
dest.write_text(json.dumps(report,indent=2),encoding="utf8")
print("STAGE980_VIDEO_SOURCE_AUDIT",json.dumps(out),flush=True)
