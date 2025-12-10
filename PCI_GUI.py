import os
import pandas as pd
import geopandas as gpd
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
from datetime import datetime

# ========= Util & Core Logic =========

def log(msg, box=None):
    print(msg)
    if box is not None:
        box.insert(tk.END, str(msg) + "\n")
        box.see(tk.END)
        box.update()

def haversine_m(lon1, lat1, lon2, lat2):
    """Haversine distance in meters. Inputs in degrees (array/pandas series ok)."""
    R = 6371008.8  # meters
    lon1 = np.radians(lon1.astype(float))
    lat1 = np.radians(lat1.astype(float))
    lon2 = np.radians(lon2.astype(float))
    lat2 = np.radians(lat2.astype(float))
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return R * c

def nearest_only(long_series, lat_series):
    """Return per-row min distance to another row (meters), ignoring self & ~0 m duplicates."""
    lon = long_series.astype(float).to_numpy()
    lat = lat_series.astype(float).to_numpy()
    n = len(lon)
    if n < 2:
        return np.full(n, np.nan, dtype="float64")

    lon_r = np.radians(lon); lat_r = np.radians(lat)
    dlon = lon_r[:, None] - lon_r[None, :]
    dlat = lat_r[:, None] - lat_r[None, :]
    a = np.sin(dlat/2.0)**2 + np.cos(lat_r[:, None]) * np.cos(lat_r[None, :]) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    dist = 6371008.8 * c

    np.fill_diagonal(dist, np.inf)
    dist[dist < 1e-6] = np.inf  # ~0 m duplicates

    nn = np.min(dist, axis=1)
    nn[np.isinf(nn)] = np.nan
    return np.round(nn, 2)

def pairwise_metrics(long_series, lat_series):
    """
    Pairwise metrics dalam satu grup:
      - min distance per-row ke tetangga terdekat (abaikan self & ~0 m)
      - index tetangga terdekat
      - min pairwise distance grup (positive)
      - max pairwise distance grup (diameter)
    """
    lon = long_series.astype(float).to_numpy()
    lat = lat_series.astype(float).to_numpy()
    n = len(lon)
    if n < 2:
        return (np.full(n, np.nan), np.full(n, -1),
                np.nan, 0.0)

    lon_r = np.radians(lon); lat_r = np.radians(lat)
    dlon = lon_r[:, None] - lon_r[None, :]
    dlat = lat_r[:, None] - lat_r[None, :]
    a = np.sin(dlat/2.0)**2 + np.cos(lat_r[:, None]) * np.cos(lat_r[None, :]) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    dist = 6371008.8 * c  # meters

    # Nearest neighbor per row
    dist_nn = dist.copy()
    np.fill_diagonal(dist_nn, np.inf)
    dist_nn[dist_nn < 1e-6] = np.inf
    nn_idx = np.argmin(dist_nn, axis=1)
    nn_dist = dist_nn[np.arange(n), nn_idx]
    nn_dist[np.isinf(nn_dist)] = np.nan

    # Group min (ignore ~0m)
    upper = np.triu(dist, k=1)
    upper_pos = np.where(upper < 1e-6, np.nan, upper)
    min_pair = np.nanmin(upper_pos) if not np.isnan(upper_pos).all() else np.nan

    # Group diameter (max pair)
    dist_no_diag = dist.copy()
    np.fill_diagonal(dist_no_diag, np.nan)
    diameter = 0.0 if np.isnan(dist_no_diag).all() else np.nanmax(dist_no_diag)

    return (np.round(nn_dist, 2), nn_idx,
            np.round(min_pair, 2) if not np.isnan(min_pair) else np.nan,
            np.round(diameter, 2))

def insert_col_at(df: pd.DataFrame, series: pd.Series, col_name: str, index: int) -> pd.DataFrame:
    """Insert/replace a column at a given index (0-based)."""
    cols = list(df.columns)
    if col_name in cols:
        cols.remove(col_name)
        df = df[cols]
    index = max(0, min(index, len(df.columns)))
    left = df.iloc[:, :index]
    right = df.iloc[:, index:]
    return pd.concat([left, series.rename(col_name), right], axis=1)

def run_processing(buffer_csv, pci_csv, buffer_radius_km, calc_group_metrics=False, logbox=None):
    """
    Pipeline:
    1) Read buffer & pci CSV
    2) Buat buffer (radius meter) -> sjoin
    3) uniq_id, uniq_count, PCI_check
    4) Q: Distance_Site_to_Buffer (buffer <-> site)
    5) R: Distance_same_uniq_m (nearest neighbor dalam uniq_id; hanya baris problem)
    6) (opsional) S,T,U,V: ringkasan grup jika calc_group_metrics=True
    7) Simpan output XLSX (timestamp) di folder yang sama dengan PCI Template
    8) Simpan GPKG (spasial)
    """
    if not buffer_csv or not os.path.exists(buffer_csv):
        raise FileNotFoundError("Buffer CSV tidak ditemukan.")
    if not pci_csv or not os.path.exists(pci_csv):
        raise FileNotFoundError("PCI Template CSV tidak ditemukan.")

    base_dir = os.path.dirname(os.path.abspath(pci_csv))
    output_dir = os.path.join(base_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    try:
        radius_km = float(buffer_radius_km)
        radius_m = radius_km * 1000.0
    except Exception:
        raise ValueError("Buffer radius (KM) tidak valid.")

    log(f"📄 Buffer CSV: {buffer_csv}", logbox)
    log(f"📄 PCI Template CSV: {pci_csv}", logbox)
    log(f"📏 Buffer radius: {radius_km} km ({radius_m:.0f} m)", logbox)
    log(f"📂 Output dir: {output_dir}", logbox)
    log(f"⚙️  Ringkasan grup S/T/U/V: {'ON' if calc_group_metrics else 'OFF'}", logbox)

    # Read CSVs
    buffer_df = pd.read_csv(buffer_csv)
    cdd_df = pd.read_csv(pci_csv)

    # Check columns
    for c in ["Longitude", "Latitude"]:
        if c not in buffer_df.columns:
            raise KeyError(f"Kolom '{c}' tidak ada di buffer CSV.")
    for c in ["Long", "Lat", "earfcn", "PCI"]:
        if c not in cdd_df.columns:
            raise KeyError(f"Kolom '{c}' tidak ada di PCI Template CSV.")

    # Site column for neighbor label
    possible_site_cols = ["Site ID", "Site_ID", "SiteID", "Site Name", "Site_Name"]
    site_col = next((c for c in possible_site_cols if c in cdd_df.columns), None)
    if site_col is None:
        site_col = "Site ID" if "Site ID" in cdd_df.columns else cdd_df.columns[0]

    # GeoDataFrames
    buffer_gdf = gpd.GeoDataFrame(
        buffer_df,
        geometry=gpd.points_from_xy(buffer_df["Longitude"], buffer_df["Latitude"]),
        crs="EPSG:4326"
    )
    cdd_gdf = gpd.GeoDataFrame(
        cdd_df,
        geometry=gpd.points_from_xy(cdd_df["Long"], cdd_df["Lat"]),
        crs="EPSG:4326"
    )

    # Buffers
    log("🛠 Membuat buffer...", logbox)
    buffer_proj = buffer_gdf.to_crs(epsg=3857)
    buffer_proj["geometry"] = buffer_proj.geometry.buffer(radius_m)
    buffer_poly = buffer_proj.to_crs(epsg=4326)

    # sjoin
    log("🔎 Spatial join (PCI ∩ buffer)...", logbox)
    intersect_gdf = gpd.sjoin(cdd_gdf, buffer_poly, how="inner", predicate="intersects")

    # uniq_id & checks
    log("🧬 uniq_id + uniq_count + PCI_check...", logbox)
    if "buffer" not in intersect_gdf.columns:
        raise KeyError("Kolom 'buffer' tidak ditemukan setelah join; pastikan ada di buffer CSV.")
    intersect_gdf["uniq_id"] = (
        intersect_gdf["buffer"].astype(str) + "_" +
        intersect_gdf["earfcn"].astype(str) + "_" +
        intersect_gdf["PCI"].astype(str)
    )
    uniq_counts = intersect_gdf["uniq_id"].value_counts().to_dict()
    intersect_gdf["uniq_count"] = intersect_gdf["uniq_id"].map(uniq_counts)
    intersect_gdf["PCI_check"] = intersect_gdf["uniq_count"].apply(lambda x: "PCI Not OK" if x > 1 else "")

    # Q: Distance_Site_to_Buffer (buffer ↔ site)
    log("📏 Hitung Distance_Site_to_Buffer (buffer ↔ site, kolom Q)...", logbox)
    need_cols = ["Long", "Lat", "Longitude", "Latitude"]
    miss = [c for c in need_cols if c not in intersect_gdf.columns]
    if miss:
        raise KeyError(f"Kolom hilang untuk Distance_Site_to_Buffer: {miss}")
    mask_ok = ~(intersect_gdf["Long"].isna() | intersect_gdf["Lat"].isna() |
                intersect_gdf["Longitude"].isna() | intersect_gdf["Latitude"].isna())
    dist_buf_site = pd.Series(np.nan, index=intersect_gdf.index, dtype="float64")
    if mask_ok.any():
        dist_buf_site.loc[mask_ok] = np.round(
            haversine_m(
                intersect_gdf.loc[mask_ok, "Long"],
                intersect_gdf.loc[mask_ok, "Lat"],
                intersect_gdf.loc[mask_ok, "Longitude"],
                intersect_gdf.loc[mask_ok, "Latitude"],
            ),
            2,
        )

    # R: nearest neighbor distance dalam uniq_id (hanya baris problem)
    log("📏 Hitung Distance_same_uniq_m (kolom R)...", logbox)
    dist_same_uniq = pd.Series(np.nan, index=intersect_gdf.index, dtype="float64")
    mask_problem = (intersect_gdf["PCI_check"].eq("PCI Not OK")) | (intersect_gdf["uniq_count"] > 1)
    if mask_problem.any():
        for uid, idx in intersect_gdf.loc[mask_problem].groupby("uniq_id").groups.items():
            sub = intersect_gdf.loc[idx, ["Long", "Lat"]]
            if sub["Long"].isna().any() or sub["Lat"].isna().any():
                continue
            dist_same_uniq.loc[idx] = nearest_only(sub["Long"], sub["Lat"])

    # === Opsional: S,T,U,V ===
    nearest_site = None
    nearest_pair_min = None
    cluster_diam = None
    pair_ok = None

    if calc_group_metrics:
        log("📊 Hitung ringkasan grup S/T/U/V...", logbox)
        nearest_site = pd.Series("", index=intersect_gdf.index, dtype="object")
        nearest_pair_min = pd.Series(np.nan, index=intersect_gdf.index, dtype="float64")
        cluster_diam = pd.Series(np.nan, index=intersect_gdf.index, dtype="float64")

        groups = intersect_gdf.loc[mask_problem].groupby("uniq_id").groups
        for uid, idx in groups.items():
            sub = intersect_gdf.loc[idx, ["Long", "Lat", site_col]]
            if sub["Long"].isna().any() or sub["Lat"].isna().any():
                continue
            nn_dist, nn_idx, g_min_pair, g_diam = pairwise_metrics(sub["Long"], sub["Lat"])
            # reuse nn_dist also to ensure consistency with R
            dist_same_uniq.loc[idx] = nn_dist

            # map nearest site ids
            vals_site = sub[site_col].tolist()
            nearest_site.loc[idx] = [("" if (k < 0 or k >= len(sub)) else str(vals_site[k])) for k in nn_idx]

            nearest_pair_min.loc[idx] = g_min_pair
            cluster_diam.loc[idx] = g_diam

        pair_ok = pd.Series(np.where((~dist_same_uniq.isna()) & (dist_same_uniq <= 2 * radius_m), True, False),
                            index=intersect_gdf.index)

    # Export XLSX (Q..V sesuai opsi) + timestamp
    log("💾 Menyimpan XLSX...", logbox)
    df_export = intersect_gdf.drop(columns="geometry").copy()
    df_export = insert_col_at(df_export, dist_buf_site, "Distance_Site_to_Buffer", 16)                       # Q
    df_export = insert_col_at(df_export, dist_same_uniq, "Distance_same_uniq_m", 17)            # R

    # only insert S..V if requested
    if calc_group_metrics:
        df_export = insert_col_at(df_export, nearest_site, "Nearest_same_uniq_SiteID", 18)      # S
        df_export = insert_col_at(df_export, nearest_pair_min, "Nearest_pair_m_per_uniq", 19)   # T
        df_export = insert_col_at(df_export, cluster_diam, "Cluster_diameter_m", 20)            # U
        df_export = insert_col_at(df_export, pair_ok, "Pair_ok_vs_buffer", 21)                  # V

    xlsx_path = os.path.join(output_dir, f"intersect_result_with_check_1_{ts}.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as w:
        df_export.to_excel(w, index=False, sheet_name="Result")

    # Spasial (tambahkan atribut yang tersedia)
    intersect_gdf["Distance_Site_to_Buffer"] = dist_buf_site
    intersect_gdf["Distance_same_uniq_m"] = dist_same_uniq
    if calc_group_metrics:
        intersect_gdf["Nearest_same_uniq_SiteID"] = nearest_site
        intersect_gdf["Nearest_pair_m_per_uniq"] = nearest_pair_min
        intersect_gdf["Cluster_diameter_m"] = cluster_diam
        intersect_gdf["Pair_ok_vs_buffer"] = pair_ok

    gpkg_path = os.path.join(output_dir, f"intersect_result_{ts}.gpkg")
    intersect_gdf.to_file(gpkg_path, layer="intersection", driver="GPKG")

    log("✅ Selesai!", logbox)
    log(f"   - {xlsx_path}", logbox)
    log(f"   - {gpkg_path}", logbox)
    return xlsx_path, gpkg_path

# ========= GUI =========

class PCIApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PCI Checking and Audit")
        self.geometry("840x610")

        # set window icon if pci.ico is present (Windows .ico)
        try:
            ico_path = os.path.join(os.path.dirname(__file__), "pci.ico")
            if os.path.exists(ico_path):
                self.iconbitmap(ico_path)
        except Exception:
            pass

        self.buffer_path_var = tk.StringVar(value="")
        self.pci_path_var = tk.StringVar(value="")
        self.radius_km_var = tk.StringVar(value="2")  # default 2 KM
        self.calc_group_var = tk.BooleanVar(value=False)  # default OFF untuk kecepatan

        # --- FOOTER (credits) ---
        style = ttk.Style(self)
        style.configure("Footer.TLabel", font=("Arial", 8), foreground="gray")
        footer = ttk.Frame(self)
        footer.pack(side="bottom", fill="x")
        credit = ttk.Label(footer, text="© 2025 - baim.muhammad@gmail.com", style="Footer.TLabel")
        credit.pack(side="right", padx=10, pady=6)

        # --- HEADER (form) ---
        header = ttk.Frame(self, padding=(12, 10, 12, 0))
        header.pack(side="top", fill="x")

        row = 0
        ttk.Label(header, text="Buffer CSV (Longitude, Latitude, buffer)").grid(row=row, column=0, sticky="w")
        ent_buf = ttk.Entry(header, textvariable=self.buffer_path_var, width=82)
        ent_buf.grid(row=row, column=1, sticky="we", padx=6, pady=(0,6))
        ttk.Button(header, text="Browse...", command=self.pick_buffer).grid(row=row, column=2, sticky="w")
        row += 1

        ttk.Label(header, text="PCI Template CSV (Long, Lat, earfcn, PCI)").grid(row=row, column=0, sticky="w")
        ent_pci = ttk.Entry(header, textvariable=self.pci_path_var, width=82)
        ent_pci.grid(row=row, column=1, sticky="we", padx=6, pady=(0,6))
        ttk.Button(header, text="Browse...", command=self.pick_pci).grid(row=row, column=2, sticky="w")
        row += 1

        ttk.Label(header, text="Buffer Radius (KM)").grid(row=row, column=0, sticky="w")
        ent_rad = ttk.Entry(header, textvariable=self.radius_km_var, width=10)
        ent_rad.grid(row=row, column=1, sticky="w", padx=6, pady=(0,2))
        row += 1

        # Checkbox untuk ringkasan S/T/U/V
        #chk = ttk.Checkbutton(header, text="Hitung kolom S,T,U,V (ringkasan grup)",
                             # variable=self.calc_group_var)
        #chk.grid(row=row, column=1, sticky="w", padx=6, pady=(0,8))
       # row += 1

        ttk.Button(header, text="Run", command=self.on_run).grid(row=row, column=0, pady=(4, 8), sticky="w")
        ttk.Button(header, text="Open Output Folder", command=self.open_output_folder).grid(row=row, column=1, pady=(4, 8), sticky="w")

        header.columnconfigure(1, weight=1)

        # --- LOG AREA (expand) ---
        logwrap = ttk.Frame(self, padding=(12, 0, 12, 0))
        logwrap.pack(side="top", fill="both", expand=True)
        self.logbox = ScrolledText(logwrap, height=18)
        self.logbox.pack(fill="both", expand=True)

        self.last_output_dir = None

    def pick_buffer(self):
        path = filedialog.askopenfilename(
            title="Pilih Buffer CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.buffer_path_var.set(path)

    def pick_pci(self):
        path = filedialog.askopenfilename(
            title="Pilih PCI Template CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.pci_path_var.set(path)

    def open_output_folder(self):
        if self.last_output_dir and os.path.isdir(self.last_output_dir):
            try:
                if os.name == 'nt':
                    os.startfile(self.last_output_dir)
                elif os.name == 'posix':
                    os.system(f'xdg-open "{self.last_output_dir}"')
                else:
                    os.system(f'open "{self.last_output_dir}"')
            except Exception as e:
                messagebox.showerror("Error", f"Tidak bisa buka folder:\n{e}")
        else:
            messagebox.showinfo("Info", "Output folder belum ada. Jalankan proses dulu.")

    def on_run(self):
        buf = self.buffer_path_var.get().strip()
        pci = self.pci_path_var.get().strip()
        rad = self.radius_km_var.get().strip()
        calc_group = self.calc_group_var.get()

        self.logbox.delete("1.0", tk.END)
        try:
            xlsx_path, gpkg_path = run_processing(buf, pci, rad, calc_group_metrics=calc_group, logbox=self.logbox)
            self.last_output_dir = os.path.dirname(xlsx_path)
            messagebox.showinfo("Selesai", f"Proses selesai.\n\nXLSX:\n{xlsx_path}\n\nGPKG:\n{gpkg_path}")
        except Exception as e:
            log(f"❌ Error: {e}", self.logbox)
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = PCIApp()
    app.mainloop()
