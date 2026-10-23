import streamlit as st
import pandas as pd
import folium
from folium import plugins
import streamlit.components.v1 as components
import os, struct, math
from collections import defaultdict, deque

st.set_page_config(page_title="한국 하천·유역 수문망 분석 시스템", layout="wide")

st.markdown("""
    <div style="background-color:#F8FAFC; padding:16px 24px; border-radius:8px; border-bottom:2px solid #E2E8F0; margin-bottom:18px;">
        <h2 style="color:#0F172A; margin:0; font-size:24px; font-weight:700;">🌊 한국 5대 수계·유역 수문망 분석 시스템</h2>
        <span style="color:#64748B; font-size:13px;">Korea River Network & Catchment Tracing System (흐름 애니메이션 및 데이터 다운로드 지원)</span>
    </div>
""", unsafe_allow_html=True)

# 1. 고정밀 좌표 변환
def tm_to_wgs84(x, y, lon_0=127.5, lat_0=38.0, k_0=0.9996, x_0=1000000.0, y_0=2000000.0):
    a = 6378137.0
    f = 1 / 298.257222101
    b = a * (1 - f)
    e2 = (a**2 - b**2) / (a**2)
    e_prime2 = (a**2 - b**2) / (b**2)
    rad_lat_0, rad_lon_0 = math.radians(lat_0), math.radians(lon_0)
    e4, e6 = e2**2, e2**3
    A0 = 1 - e2/4 - 3*e4/64 - 5*e6/256
    A2 = 3/8 * (e2 + e4/4 + 15*e6/128)
    A4 = 15/256 * (e4 + 3*e6/4)
    A6 = 35*e6/3072
    M0 = a * (A0*rad_lat_0 - A2*math.sin(2*rad_lat_0) + A4*math.sin(4*rad_lat_0) - A6*math.sin(6*rad_lat_0))
    M = M0 + (y - y_0) / k_0
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))
    mu = M / (a * A0)
    J1, J2, J3, J4 = 3*e1/2 - 27*e1**3/32, 21*e1**2/16 - 55*e1**4/32, 151*e1**3/96, 1097*e1**4/512
    fp = mu + J1*math.sin(2*mu) + J2*math.sin(4*mu) + J3*math.sin(6*mu) + J4*math.sin(8*mu)
    sin_fp, cos_fp, tan_fp = math.sin(fp), math.cos(fp), math.tan(fp)
    C1 = e_prime2 * (cos_fp**2)
    T1 = tan_fp**2
    N1 = a / math.sqrt(1 - e2 * sin_fp**2)
    R1 = a * (1 - e2) / math.pow(1 - e2 * sin_fp**2, 1.5)
    D = (x - x_0) / (N1 * k_0)
    lat = fp - (N1 * tan_fp / R1) * (D**2/2 - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*e_prime2)*D**4/24 + (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*e_prime2 - 3*C1**2)*D**6/720)
    lon = rad_lon_0 + (D - (1 + 2*T1 + C1)*D**3/6 + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*e_prime2 + 24*T1**2)*D**5/120) / cos_fp
    return math.degrees(lat), math.degrees(lon)

def web_mercator_to_wgs84(x, y):
    lon = math.degrees(x / 6378137.0)
    lat = math.degrees(2 * math.atan(math.exp(y / 6378137.0)) - math.pi / 2)
    return lat, lon

# 2. 데이터 로드 및 파싱 함수 모음
@st.cache_data
def load_krf_data():
    for p in ["data/표준유역 네트워크_KRF수정.xlsx", "표준유역 네트워크_KRF수정.xlsx"]:
        if os.path.exists(p):
            df = pd.read_excel(p, sheet_name="Sheet1")
            up_graph = defaultdict(list)
            for _, r in df.iterrows():
                c = str(r['표준유역코드']).strip()
                d = str(r['직접하류_표준유역코드']).strip()
                if d and d not in ['nan', 'OUTLET']:
                    up_graph[d].append(c)
            return df, up_graph
    return None, None

@st.cache_data
def load_dams():
    shp_path = "data/dm_obsvtr.shp" if os.path.exists("data/dm_obsvtr.shp") else "dm_obsvtr.shp"
    dbf_path = shp_path.replace(".shp", ".dbf")
    dams = []
    if os.path.exists(shp_path) and os.path.exists(dbf_path):
        with open(dbf_path, "rb") as f_dbf:
            header = f_dbf.read(32)
            n_rec, h_len, r_len = struct.unpack('<IHH', header[4:12])
            f_dbf.seek(h_len)
            names = [f_dbf.read(r_len)[255:509].decode("utf-8", errors="ignore").strip() for _ in range(n_rec)]
        with open(shp_path, "rb") as f_shp:
            f_shp.seek(100)
            for i in range(len(names)):
                f_shp.read(8)
                f_shp.read(4)
                x, y = struct.unpack("<dd", f_shp.read(16))
                lat, lon = web_mercator_to_wgs84(x, y)
                dams.append({"name": names[i], "lat": round(lat, 5), "lon": round(lon, 5)})
    return dams

def get_river_polygons(target_name):
    paths_to_check = [
        ("data/00.하천망도_국가.shp", "data/00.하천망도_국가.dbf"),
        ("00.하천망도_국가.shp", "00.하천망도_국가.dbf"),
        ("data/00.하천망도_지방.shp", "data/00.하천망도_지방.dbf"),
        ("00.하천망도_지방.shp", "00.하천망도_지방.dbf"),
    ]
    found_polygons = []
    for shp_p, dbf_p in paths_to_check:
        shx_p = shp_p.replace(".shp", ".shx")
        if not (os.path.exists(shp_p) and os.path.exists(dbf_p) and os.path.exists(shx_p)):
            continue
        matches = []
        with open(dbf_p, "rb") as f_dbf:
            header = f_dbf.read(32)
            n_rec, h_len, r_len = struct.unpack('<IHH', header[4:12])
            f_dbf.seek(h_len)
            for i in range(n_rec):
                raw = f_dbf.read(r_len)
                r_name = raw[21:275].decode("utf-8", errors="ignore").strip()
                if target_name == r_name or target_name in r_name:
                    matches.append((i, r_name))
        if not matches:
            continue
        with open(shx_p, "rb") as f_shx, open(shp_p, "rb") as f_shp:
            for idx, r_name in matches:
                f_shx.seek(100 + idx * 8)
                offset, _ = struct.unpack(">II", f_shx.read(8))
                f_shp.seek(offset * 2 + 8)
                f_shp.read(4)
                f_shp.read(32)
                num_parts, num_points = struct.unpack("<II", f_shp.read(8))
                parts = struct.unpack(f"<{num_parts}I", f_shp.read(4 * num_parts))
                part_bounds = list(parts) + [num_points]
                raw_pts = [struct.unpack("<dd", f_shp.read(16)) for _ in range(num_points)]
                for p in range(num_parts):
                    p_pts = raw_pts[part_bounds[p]:part_bounds[p+1]]
                    step = max(1, len(p_pts) // 600)
                    coords = []
                    for j in range(0, len(p_pts), step):
                        lat, lon = tm_to_wgs84(p_pts[j][0], p_pts[j][1])
                        coords.append([round(lat, 6), round(lon, 6)])
                    if coords:
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                        found_polygons.append({"name": r_name, "coords": coords})
        if found_polygons:
            break
    return found_polygons

@st.cache_data
def get_basin_polygons_from_shp(target_prefix):
    shp_paths = ["data/mgnl_sbsn_bsn.shp", "mgnl_sbsn_bsn.shp"]
    target_shp = next((p for p in shp_paths if os.path.exists(p) and os.path.exists(p.replace(".shp", ".dbf"))), None)
    if not target_shp: return None
        
    dbf_p = target_shp.replace(".shp", ".dbf")
    shx_p = target_shp.replace(".shp", ".shx")
    prefix_str = str(target_prefix).strip()
    
    matches = []
    with open(dbf_p, "rb") as f_dbf:
        header = f_dbf.read(32)
        n_rec, h_len, r_len = struct.unpack('<IHH', header[4:12])
        f_dbf.seek(h_len)
        for i in range(n_rec):
            raw = f_dbf.read(r_len)
            c = raw[1:255].decode("utf-8", errors="ignore").strip()
            name = raw[255:509].decode("utf-8", errors="ignore").strip()
            if len(prefix_str) == 6 and c == prefix_str:
                matches.append((i, c, name))
                break
            elif len(prefix_str) <= 4 and c.startswith(prefix_str):
                matches.append((i, c, name))
                
    if not matches: return None
        
    polys = []
    with open(shx_p, "rb") as f_shx, open(target_shp, "rb") as f_shp:
        for idx, c, name in matches:
            f_shx.seek(100 + idx * 8)
            offset, _ = struct.unpack(">II", f_shx.read(8))
            f_shp.seek(offset * 2 + 8 + 4)
            xmin, ymin, xmax, ymax = struct.unpack("<dddd", f_shp.read(32))
            num_parts, num_points = struct.unpack("<II", f_shp.read(8))
            parts = struct.unpack(f"<{num_parts}I", f_shp.read(4 * num_parts))
            part_bounds = list(parts) + [num_points]
            raw_pts = [struct.unpack("<dd", f_shp.read(16)) for _ in range(num_points)]
            
            step = max(1, num_points // (120 if len(prefix_str) <= 4 else 300))
            center_pt = [round((ymin + ymax) / 2, 5), round((xmin + xmax) / 2, 5)]
            
            for p in range(num_parts):
                p_pts = raw_pts[part_bounds[p]:part_bounds[p+1]]
                coords = [[round(p_pts[j][1], 5), round(p_pts[j][0], 5)] for j in range(0, len(p_pts), step)]
                if coords and coords[0] != coords[-1]:
                    coords.append(coords[0])
                polys.append({"code": c, "name": name, "coords": coords, "center": center_pt})
    return polys

@st.cache_data
def load_all_basin_centers():
    shp_paths = ["data/mgnl_sbsn_bsn.shp", "mgnl_sbsn_bsn.shp"]
    target_shp = next((p for p in shp_paths if os.path.exists(p) and os.path.exists(p.replace(".shp", ".dbf"))), None)
    if not target_shp: return {}

    dbf_p = target_shp.replace(".shp", ".dbf")
    shx_p = target_shp.replace(".shp", ".shx")
    
    centers = {}
    with open(dbf_p, "rb") as f_dbf:
        header = f_dbf.read(32)
        n_rec, h_len, r_len = struct.unpack('<IHH', header[4:12])
        f_dbf.seek(h_len)
        for i in range(n_rec):
            raw = f_dbf.read(r_len)
            c = raw[1:255].decode("utf-8", errors="ignore").strip()
            centers[c] = i
            
    res = {}
    with open(shx_p, "rb") as f_shx, open(target_shp, "rb") as f_shp:
        for c, idx in centers.items():
            f_shx.seek(100 + idx * 8)
            offset, _ = struct.unpack(">II", f_shx.read(8))
            f_shp.seek(offset * 2 + 8 + 4)
            xmin, ymin, xmax, ymax = struct.unpack("<dddd", f_shp.read(32))
            res[c] = [round((ymin + ymax) / 2, 5), round((xmin + xmax) / 2, 5)]
    return res

# 데이터프레임 CSV 변환 함수
@st.cache_data
def convert_df_to_csv(df):
    return df.to_csv(index=False, encoding='utf-8-sig')

df_krf, up_graph = load_krf_data()
all_dams = load_dams()
all_centers = load_all_basin_centers()

# 5대 대권역 메타데이터
MAJOR_BASINS = {
    "한강권역": {"codes": [10, 11, 12, 13], "center": [37.60, 127.50], "zoom": 8, "mains": ["한강", "안성천", "북한강", "소양강", "임진강"]},
    "낙동강권역": {"codes": [20, 21, 22, 23, 24, 25], "center": [35.80, 128.60], "zoom": 8, "mains": ["낙동강", "금호강", "남강", "형산강", "태화강"]},
    "금강권역": {"codes": [30, 31, 32, 33], "center": [36.35, 127.10], "zoom": 8, "mains": ["금강", "갑천", "미호천", "삽교천", "만경강", "동진강"]},
    "섬진강권역": {"codes": [40, 41], "center": [35.25, 127.35], "zoom": 9, "mains": ["섬진강", "보성강", "요천"]},
    "영산강권역": {"codes": [50, 51, 52, 53], "center": [34.95, 126.70], "zoom": 9, "mains": ["영산강", "황룡강", "광주천", "탐진강"]},
    "제주도권역": {"codes": [60], "center": [33.38, 126.55], "zoom": 10, "mains": []}
}

# 사이드바 UI 구성
st.sidebar.markdown("### 🌐 1. 권역 선택 (5대 대권역)")
sel_major = st.sidebar.selectbox("1) 대권역 선택", list(MAJOR_BASINS.keys()), index=0)
major_meta = MAJOR_BASINS[sel_major]
sub_codes = major_meta["codes"]
df_major = df_krf[df_krf['대권역코드'].isin(sub_codes)]

available_jung = sorted(df_major['중권역코드'].unique())
jung_options = ["🔍 [대권역 전체 보기]"]
for jc in available_jung:
    grp = df_major[df_major['중권역코드'] == jc]
    main_r = grp['주요하천'].dropna().tolist()
    r_label = main_r[0].split(';')[0].strip() if main_r else grp.iloc[0]['표준유역명']
    jung_options.append(f"{jc} ({r_label}권 / {grp.iloc[0]['표준유역명']}~{grp.iloc[-1]['표준유역명']})")

sel_jung_label = st.sidebar.selectbox("2) 중권역 선택", jung_options, index=0)

view_level = "대권역"
sel_jung_code = None
selected_basin_code = None

if sel_jung_label != "🔍 [대권역 전체 보기]":
    sel_jung_code = int(sel_jung_label.split()[0])
    df_jung = df_major[df_major['중권역코드'] == sel_jung_code]
    
    basin_options = ["🔍 [중권역 전체 보기]"] + [f"{r['표준유역코드']} - {r['표준유역명']}" for _, r in df_jung.iterrows()]
    sel_basin_label = st.sidebar.selectbox("3) 표준유역 선택", basin_options, index=0)
    
    if sel_basin_label == "🔍 [중권역 전체 보기]":
        view_level = "중권역"
    else:
        view_level = "표준유역"
        selected_basin_code = sel_basin_label.split(" - ")[0].strip()

st.sidebar.markdown("---")
st.sidebar.markdown("### 🗺️ 2. 지도 표시 옵션")
show_dams = st.sidebar.checkbox("🏢 댐 (위치 마커)", value=True)
show_river_net = st.sidebar.checkbox("🌊 하천 네트워크", value=True)
show_upstream_basin = st.sidebar.checkbox("⛰️ 상류유역 (집수구역 면 및 마크)", value=True)
show_related_basin = st.sidebar.checkbox("🗺️ 관련 표준유역 (현재 영역)", value=True)
show_flow_dir = st.sidebar.checkbox("↪️ 각 표준유역 물 흐름 방향 (네트워크)", value=True)

# KRF 분석 데이터 계산 및 다운로드용 데이터프레임 조립
DAM_CONT = {"100319": 100320, "101207": 101208}
full_names, full_codes = [], []
upstream_codes, upstream_names = [], []

center = major_meta["center"]
zoom_val = major_meta["zoom"]
display_title = ""

# 다운로드용 DataFrame 초기화
df_down_export = pd.DataFrame()
df_list_export = pd.DataFrame()
list_export_filename = ""
list_button_label = ""

if view_level == "대권역":
    display_title = f"{sel_major} 전체 현황"
    dae_summary = []
    for jc in available_jung:
        g = df_major[df_major['중권역코드'] == jc]
        mr = g['주요하천'].dropna().tolist()
        rep_r = mr[0].split(';')[0].strip() if mr else g.iloc[0]['표준유역명']
        dae_summary.append({
            "중권역코드": jc,
            "대표 하천/수계": rep_r,
            "포함 표준유역 수": len(g),
            "기점 ~ 종점 유역": f"{g.iloc[0]['표준유역명']} ➔ {g.iloc[-1]['표준유역명']}"
        })
    df_list_export = pd.DataFrame(dae_summary)
    list_export_filename = f"중권역목록_{sel_major}.csv"
    list_button_label = "⬇️ 소속 중권역 목록 다운로드"

elif view_level == "중권역":
    display_title = f"중권역 {sel_jung_label.split()[0]} (총 {len(df_jung)}개 유역)"
    all_j_codes = set(df_jung['표준유역코드'].astype(str))
    outlet_basin = df_jung.iloc[-1]
    for _, r in df_jung.iterrows():
        if str(r['직접하류_표준유역코드']).strip() not in all_j_codes:
            outlet_basin = r
            break
            
    p_names = [p.strip() for p in str(outlet_basin.get("하류전체경로_명칭", "")).split("→")]
    p_codes = [c.strip() for c in str(outlet_basin.get("하류전체경로_코드", "")).split("→")]
    if len(p_codes) >= 2 and p_codes[-2] in DAM_CONT:
        next_c = DAM_CONT[p_codes[-2]]
        dm_row = df_krf[df_krf["표준유역코드"] == next_c]
        if len(dm_row) > 0:
            p_names = p_names[:-1] + [p.strip() for p in str(dm_row.iloc[0].get("하류전체경로_명칭", "")).split("→")]
            p_codes = p_codes[:-1] + [c.strip() for c in str(dm_row.iloc[0].get("하류전체경로_코드", "")).split("→")]
    full_names, full_codes = p_names, p_codes

    df_list_export = df_jung[['표준유역코드', '표준유역명', '국가하천', '지방하천', '직접하류_표준유역명']].copy()
    list_export_filename = f"표준유역목록_{sel_jung_code}.csv"
    list_button_label = "⬇️ 소속 표준유역 목록 다운로드"

elif view_level == "표준유역":
    m_row = df_krf[df_krf["표준유역코드"].astype(str) == selected_basin_code].iloc[0]
    display_title = f"{m_row['표준유역명']} ({selected_basin_code})"
    
    p_names = [p.strip() for p in str(m_row.get("하류전체경로_명칭", "")).split("→")]
    p_codes = [c.strip() for c in str(m_row.get("하류전체경로_코드", "")).split("→")]
    if len(p_codes) >= 2 and p_codes[-2] in DAM_CONT:
        next_c = DAM_CONT[p_codes[-2]]
        dm_row = df_krf[df_krf["표준유역코드"] == next_c]
        if len(dm_row) > 0:
            p_names = p_names[:-1] + [p.strip() for p in str(dm_row.iloc[0].get("하류전체경로_명칭", "")).split("→")]
            p_codes = p_codes[:-1] + [c.strip() for c in str(dm_row.iloc[0].get("하류전체경로_코드", "")).split("→")]
    full_names, full_codes = p_names, p_codes

    if up_graph is not None:
        visited = set()
        queue = deque([selected_basin_code])
        while queue:
            curr = queue.popleft()
            for up in up_graph.get(curr, []):
                if up not in visited:
                    visited.add(up)
                    queue.append(up)
        upstream_codes = sorted(list(visited))
        up_df = df_krf[df_krf['표준유역코드'].astype(str).isin(upstream_codes)]
        upstream_names = up_df['표준유역명'].tolist()
        
        if upstream_codes:
            up_records = [{"순번": i+1, "유역코드": c, "유역명": df_krf[df_krf['표준유역코드'].astype(str)==c].iloc[0]['표준유역명']} for i, c in enumerate(upstream_codes)]
            df_list_export = pd.DataFrame(up_records)
            
    list_export_filename = f"상류기여유역_{selected_basin_code}.csv"
    list_button_label = "⬇️ 상류 기여구역 다운로드"

# 공통 하류경로 다운로드 데이터프레임
if full_names:
    df_down_export = pd.DataFrame([{"단계": i+1, "표준유역코드": full_codes[i], "표준유역명": full_names[i]} for i in range(len(full_names))])

# 메인 화면 레이아웃
col_main, col_sub = st.columns([2.1, 1.0], gap="large")

with col_main:
    # 1. 지도 중심 및 줌 계산
    if view_level == "표준유역":
        local_polys = get_basin_polygons_from_shp(selected_basin_code)
        basin_name = m_row['표준유역명']
        primary_river = str(m_row.get('주요하천', '')).split(';')[0].strip()
        river_polys = get_river_polygons(primary_river) or get_river_polygons(basin_name)
        
        if local_polys:
            all_lats = [pt[0] for poly in local_polys for pt in poly['coords']]
            all_lons = [pt[1] for poly in local_polys for pt in poly['coords']]
            center = [(min(all_lats)+max(all_lats))/2, (min(all_lons)+max(all_lons))/2]
            zoom_val = 11
        elif river_polys:
            all_lats = [pt[0] for p in river_polys for pt in p["coords"]]
            all_lons = [pt[1] for p in river_polys for pt in p["coords"]]
            center = [(min(all_lats)+max(all_lats))/2, (min(all_lons)+max(all_lons))/2]
            zoom_val = 11
            
    elif view_level == "중권역":
        jung_polys = get_basin_polygons_from_shp(sel_jung_code)
        if jung_polys:
            all_lats = [pt[0] for poly in jung_polys for pt in poly['coords']]
            all_lons = [pt[1] for poly in jung_polys for pt in poly['coords']]
            center = [(min(all_lats)+max(all_lats))/2, (min(all_lons)+max(all_lons))/2]
            zoom_val = 10
        else:
            center = major_meta["center"]
            zoom_val = 10

    m = folium.Map(location=center, zoom_start=zoom_val, tiles="OpenStreetMap")

    # 옵션 A: 상류유역(집수구역)
    if show_upstream_basin and view_level == "표준유역" and upstream_codes:
        added_centers = set()
        for u_c in upstream_codes:
            u_polys = get_basin_polygons_from_shp(u_c)
            if u_polys:
                for poly in u_polys:
                    folium.Polygon(
                        locations=poly['coords'], color="#059669", weight=1.2, fill=True,
                        fill_color="#10B981", fill_opacity=0.25, tooltip=f"상류 기여유역: {poly['name']} ({poly['code']})"
                    ).add_to(m)
                    c_lat, c_lon = poly["center"]
                    pos_key = f"{c_lat:.3f}_{c_lon:.3f}"
                    if pos_key not in added_centers:
                        added_centers.add(pos_key)
                        folium.Marker(
                            [c_lat, c_lon], popup=f"⛰️ [상류 집수구역] {poly['name']}", tooltip=f"상류: {poly['name']}",
                            icon=folium.Icon(color="green", icon="tree", prefix="fa")
                        ).add_to(m)

    # 옵션 B: 관련 표준유역
    if show_related_basin:
        if view_level == "표준유역" and local_polys:
            for poly in local_polys:
                folium.Polygon(
                    locations=poly['coords'], color="#DC2626", weight=2.5, fill=True,
                    fill_color="#EF4444", fill_opacity=0.35, tooltip=f"선택 표준유역: {poly['name']} ({poly['code']})"
                ).add_to(m)
        elif view_level == "중권역" and 'jung_polys' in locals() and jung_polys:
            for poly in jung_polys:
                folium.Polygon(
                    locations=poly['coords'], color="#EA580C", weight=1.5, fill=True,
                    fill_color="#F97316", fill_opacity=0.25, tooltip=f"중권역 표준유역: {poly['name']} ({poly['code']})"
                ).add_to(m)

    # 옵션 C: 하천 네트워크
    if show_river_net:
        if view_level == "대권역":
            for r_name in major_meta["mains"]:
                polys = get_river_polygons(r_name)
                for poly in polys:
                    folium.Polygon(
                        locations=poly["coords"], color="#2563EB", weight=2.0, fill=True,
                        fill_color="#60A5FA", fill_opacity=0.6, tooltip=f"권역 본류: {poly['name']}"
                    ).add_to(m)
        elif view_level == "중권역":
            for _, r in df_jung.iterrows():
                r_main = str(r['주요하천']).split(';')[0].strip() if pd.notna(r['주요하천']) else r['표준유역명']
                polys = get_river_polygons(r_main)
                for poly in polys:
                    folium.Polygon(
                        locations=poly["coords"], color="#0284C7", weight=2.0, fill=True,
                        fill_color="#38BDF8", fill_opacity=0.6, tooltip=f"하천망: {poly['name']}"
                    ).add_to(m)
        elif view_level == "표준유역":
            if river_polys:
                for poly in river_polys:
                    folium.Polygon(
                        locations=poly["coords"], color="#B91C1C", weight=2.5, fill=True,
                        fill_color="#EF4444", fill_opacity=0.8, tooltip=f"하천구역: {poly['name']}"
                    ).add_to(m)

    # 옵션 D: 댐 마커
    if show_dams and all_dams:
        for d in all_dams:
            is_relevant = any(d["name"].replace("댐", "") in nm for nm in full_names) or view_level in ["대권역", "중권역"]
            if is_relevant:
                folium.Marker(
                    [d["lat"], d["lon"]], popup=f"🏢 [댐] {d['name']} (수계 조절점)", tooltip=d['name'],
                    icon=folium.Icon(color="darkblue", icon="tint")
                ).add_to(m)

    # 옵션 E: 각 표준유역 물 흐름 방향 (네트워크)
    if show_flow_dir and all_centers:
        flow_edges = set()
        for k, v in DAM_CONT.items():
            flow_edges.add((str(k), str(v)))

        if view_level == "대권역":
            for _, r in df_major.iterrows():
                c1, c2 = str(r['표준유역코드']).strip(), str(r['직접하류_표준유역코드']).strip()
                if c2 and c2 != 'nan' and c2 != 'OUTLET':
                    flow_edges.add((c1, c2))
        elif view_level == "중권역":
            for _, r in df_jung.iterrows():
                c1, c2 = str(r['표준유역코드']).strip(), str(r['직접하류_표준유역코드']).strip()
                if c2 and c2 != 'nan' and c2 != 'OUTLET':
                    flow_edges.add((c1, c2))
            for i in range(len(full_codes)-1):
                flow_edges.add((full_codes[i], full_codes[i+1]))
        elif view_level == "표준유역":
            if upstream_codes:
                for u_c in upstream_codes:
                    u_row = df_krf[df_krf['표준유역코드'].astype(str) == u_c]
                    if len(u_row) > 0:
                        c2 = str(u_row.iloc[0]['직접하류_표준유역코드']).strip()
                        if c2 and c2 != 'nan' and c2 != 'OUTLET':
                            flow_edges.add((u_c, c2))
            for i in range(len(full_codes)-1):
                flow_edges.add((full_codes[i], full_codes[i+1]))

        flow_segments = []
        for c1, c2 in flow_edges:
            if c1 in all_centers and c2 in all_centers:
                flow_segments.append([all_centers[c1], all_centers[c2]])
        
        if flow_segments:
            plugins.AntPath(
                locations=flow_segments,
                dash_array=[10, 20],
                delay=800,
                color='#1D4ED8',
                pulse_color='#FFFFFF',
                weight=3.5,
                tooltip="🌊 유출 흐름 방향"
            ).add_to(m)

    components.html(m._repr_html_(), height=520)

    # 하단 탭 영역 (테이블 출력)
    if view_level == "대권역":
        st.markdown(f"#### 📊 {sel_major} 소속 중권역 총괄 ({len(df_list_export)}개 중권역)")
        st.dataframe(df_list_export, use_container_width=True, hide_index=True)

    elif view_level == "중권역":
        t1, t2 = st.tabs([f"🌊 1. 중권역 출구 ➔ 최종 방류 경로", f"📍 2. 중권역 내 표준유역 목록 ({len(df_list_export)}개)"])
        with t1:
            st.dataframe(df_down_export, use_container_width=True, hide_index=True)
        with t2:
            st.dataframe(df_list_export, use_container_width=True, hide_index=True)

    elif view_level == "표준유역":
        t1, t2, t3 = st.tabs([f"🌊 1. 하류 방류 경로 ({len(df_down_export)}단계)", f"⛰️ 2. 전체 상류 기여구역 ({len(df_list_export)}개 유역)", f"📍 3. 집수구역 상세"])
        with t1:
            st.dataframe(df_down_export, use_container_width=True, hide_index=True)
        with t2:
            if not df_list_export.empty:
                st.dataframe(df_list_export, use_container_width=True, hide_index=True)
            else:
                st.info("이 유역은 수계 최상류 발원지입니다.")
        with t3:
            st.dataframe(pd.DataFrame([
                {"항목": "표준유역", "내용": f"{m_row['표준유역명']} ({selected_basin_code})"},
                {"항목": "중권역", "내용": str(m_row['중권역코드'])},
                {"항목": "국가하천", "내용": str(m_row.get('국가하천', '-'))},
                {"항목": "지방하천", "내용": str(m_row.get('지방하천', '-'))},
                {"항목": "소하천", "내용": str(m_row.get('소하천', '-'))}
            ]), use_container_width=True, hide_index=True)

with col_sub:
    st.markdown("#### 📍 수계 제원 요약")
    st.markdown(f"""
        <div style="background-color:#F1F5F9; padding:14px; border-radius:8px; border:1px solid #CBD5E1; font-size:13px; line-height:1.7;">
            <b>현재 분석 레벨:</b> <span style="color:#2563EB; font-weight:700;">{view_level}</span><br>
            <b>대권역:</b> {sel_major}<br>
            <b>선택 대상:</b> {display_title}<br>
        </div>
    """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    if not df_down_export.empty:
        st.markdown(f"#### 🧭 유출 플로우 단계 ({len(full_names)}개 구간)")
        flow_box = """<div style="background-color:#FFFFFF; padding:12px; border-radius:8px; border:1px solid #E2E8F0; font-size:13px; max-height:220px; overflow-y:auto;">"""
        for idx, nm in enumerate(full_names):
            bg = "#EFF6FF" if idx == 0 else ("#E0F2FE" if idx == len(full_names)-1 else "#F8FAFC")
            col = "#1E40AF" if idx == 0 else ("#0369A1" if idx == len(full_names)-1 else "#334155")
            flow_box += f"""<div style="padding:6px 10px; background:{bg}; border-radius:6px; margin-bottom:4px; font-weight:600; color:{col};">{idx+1}. {nm}</div>"""
            if idx < len(full_names) - 1:
                flow_box += """<div style="text-align:center; color:#94A3B8; font-size:11px; margin:-2px 0;">↓</div>"""
        flow_box += "</div>"
        st.markdown(flow_box, unsafe_allow_html=True)
    else:
        st.info("권역 내 중권역 목록이 하단 테이블에 요약 표시 중입니다.")
        
    # --- 데이터 다운로드 버튼 구역 ---
    st.markdown("<br>#### 💾 분석 데이터 다운로드", unsafe_allow_html=True)
    
    if not df_down_export.empty:
        st.download_button(
            label="⬇️ 하류 방류 경로 다운로드 (CSV)",
            data=convert_df_to_csv(df_down_export),
            file_name=f"하류경로_{display_title.replace(' ', '_')}.csv",
            mime="text/csv",
            use_container_width=True
        )
        
    if not df_list_export.empty:
        st.download_button(
            label=list_button_label,
            data=convert_df_to_csv(df_list_export),
            file_name=list_export_filename,
            mime="text/csv",
            use_container_width=True
        )