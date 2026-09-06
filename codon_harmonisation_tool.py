import csv
import io
import re
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP_TITLE = "Rare Codon Analysis Tool"

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

SYNONYMS = {}
for codon, aa in GENETIC_CODE.items():
    SYNONYMS.setdefault(aa, []).append(codon)
for aa in SYNONYMS:
    SYNONYMS[aa] = sorted(SYNONYMS[aa])

EXAMPLE_SOURCE = """CGT 4.5
CGC 10.2
CGA 1.1
CGG 0.9
AGA 2.2
AGG 0.8
AAA 24.4
AAG 12.1
ATG 22.0
TTT 17.6
TTC 20.4
"""
EXAMPLE_TARGET = """CGT 8.1
CGC 20.2
CGA 3.0
CGG 2.0
AGA 1.0
AGG 0.6
AAA 10.4
AAG 29.8
ATG 23.0
TTT 10.0
TTC 30.0
"""

HELP_TEXT = """
Unterstützte Formate für Codonfrequenz-Tabellen:
1) Freitext / Kazuza-ähnlich:
   AAA 24.4
   AAG 12.1
   ...
2) Kazuza-Block mit Klammerzahlen:
   UUU 26.1(170666)  UCU 23.5(153557) ...
3) CSV/TSV mit mindestens:
   codon, frequency
   optional zusätzlich: amino_acid

Das Tool verwendet den Standard-Genetischen Code.
Score eines Codons = Frequenz(Codon) / Summe der synonymen Codonfrequenzen.
RFN (Relative Frequency Number) = Rang innerhalb einer Aminosäure,
sortiert nach Score aufsteigend (1 = seltenstes synonymes Codon).

Ziel-Matching:
- Für jede Position wird ein Target-Codon gewählt, das zur gewünschten Aminosäure passt
  und dessen Target-Score dem Source-Score möglichst ähnlich ist.
- Falls Aminosäuresequenz und Source-CDS nicht zusammenpassen, wird die Aminosäuresequenz
  als Soll-Sequenz verwendet und das Codon entsprechend korrigiert.
"""

def clean_dna(seq: str) -> str:
    seq = re.sub(r"[^ACGTUacgtu]", "", seq or "")
    return seq.upper().replace("U", "T")

def clean_aa(seq: str) -> str:
    seq = re.sub(r"[^A-Za-z\*]", "", seq or "")
    return seq.upper()

def translate_dna(seq: str) -> str:
    if len(seq) % 3 != 0:
        raise ValueError("Die Coding-Sequenz hat keine Länge, die durch 3 teilbar ist.")
    aas = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon not in GENETIC_CODE:
            raise ValueError(f"Ungültiges Codon gefunden: {codon}")
        aas.append(GENETIC_CODE[codon])
    return "".join(aas)

def parse_frequency_table(text: str) -> dict:
    text = text.strip()
    if not text:
        raise ValueError("Leere Codonfrequenz-Tabelle.")

    if "," in text or "\t" in text or ";" in text:
        parsed = parse_delimited_table(text)
        if parsed:
            return parsed

    matches = re.findall(
        r"\b([ACGTUacgtu]{3})\b\s*[:=,;\t ]+\s*([0-9]+(?:[.,][0-9]+)?)(?:\s*\(\s*[\d\s]+\))?",
        text
    )
    if not matches:
        raise ValueError("Konnte keine Codon/Frequenz-Paare erkennen.")
    freq = {}
    for codon, value in matches:
        codon = codon.upper().replace("U", "T")
        if codon in GENETIC_CODE:
            freq[codon] = float(value.replace(",", "."))
    if not freq:
        raise ValueError("Keine gültigen Codons in der Tabelle erkannt.")
    return fill_missing_codons(freq)

def parse_delimited_table(text: str) -> dict:
    for delimiter in [",", "\t", ";"]:
        try:
            reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
            if not reader.fieldnames:
                continue
            field_map = {name.strip().lower(): name for name in reader.fieldnames}
            codon_key = None
            freq_key = None
            for candidate in ["codon", "triplet", "codons"]:
                if candidate in field_map:
                    codon_key = field_map[candidate]
                    break
            for candidate in ["frequency", "freq", "per_thousand", "value", "count"]:
                if candidate in field_map:
                    freq_key = field_map[candidate]
                    break
            if not codon_key or not freq_key:
                continue

            freq = {}
            for row in reader:
                codon = (row.get(codon_key) or "").strip().upper().replace("U", "T")
                value = (row.get(freq_key) or "").strip().replace(",", ".")
                if codon in GENETIC_CODE and value:
                    freq[codon] = float(value)
            if freq:
                return fill_missing_codons(freq)
        except Exception:
            continue
    return {}

def fill_missing_codons(freq: dict) -> dict:
    full = {c: 0.0 for c in GENETIC_CODE.keys()}
    full.update(freq)
    return full

def compute_metrics(freq_table: dict) -> dict:
    result = {}
    for aa, codons in SYNONYMS.items():
        total = sum(freq_table.get(c, 0.0) for c in codons)
        sortable = []
        for codon in codons:
            value = freq_table.get(codon, 0.0)
            score = (value / total) if total > 0 else 0.0
            sortable.append((codon, value, score))
        sortable.sort(key=lambda x: (x[2], x[1], x[0]))
        for rank, (codon, value, score) in enumerate(sortable, start=1):
            result[codon] = {
                "aa": aa,
                "frequency": value,
                "score": score,
                "rfn": rank,
                "synonymous_total": total,
            }
    return result

def format_float(x):
    return f"{x:.6f}".rstrip("0").rstrip(".")

def choose_best_target_codon(desired_aa: str, source_score: float, target_metrics: dict) -> tuple[str, float]:
    candidates = SYNONYMS.get(desired_aa, [])
    if not candidates:
        raise ValueError(f"Keine Target-Codons für Aminosäure {desired_aa} gefunden.")
    ranked = sorted(
        candidates,
        key=lambda c: (
            abs(target_metrics[c]["score"] - source_score),
            -target_metrics[c]["frequency"],
            c,
        ),
    )
    best = ranked[0]
    return best, abs(target_metrics[best]["score"] - source_score)

def analyze_sequence(
    source_cds: str,
    desired_aa_seq: str,
    source_metrics: dict,
    target_metrics: dict,
):
    translated = translate_dna(source_cds)
    rows = []
    target_codons = []
    target_aas = []

    for idx, pos in enumerate(range(0, len(source_cds), 3), start=1):
        codon = source_cds[pos:pos+3]
        source_aa = GENETIC_CODE[codon]
        desired_aa = desired_aa_seq[idx - 1]
        src = source_metrics[codon]
        source_score = src["score"]

        mapped, score_delta = choose_best_target_codon(
            desired_aa=desired_aa,
            source_score=source_score,
            target_metrics=target_metrics,
        )

        target_codons.append(mapped)
        target_aas.append(desired_aa)

        notes = []
        if source_aa != desired_aa:
            notes.append(f"AA-Mismatch: Source-Codon kodiert {source_aa}, Soll-AA ist {desired_aa}.")
        notes.append(f"Target-Codon über minimalen Score-Abstand gewählt (Δ={format_float(score_delta)}).")

        rows.append(
            {
                "Position": idx,
                "AA": desired_aa,
                "Source Codon": codon,
                "Source Frequency": src["frequency"],
                "Source Score": src["score"],
                "Source RFN": src["rfn"],
                "Target Codon": mapped,
                "Target Frequency": target_metrics[mapped]["frequency"],
                "Target Score": target_metrics[mapped]["score"],
                "Target RFN": target_metrics[mapped]["rfn"],
                "Note": " ".join(notes),
            }
        )

    return rows, "".join(target_codons), "".join(target_aas), translated

class RareCodonApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1420x900")
        self.minsize(1180, 760)

        self.source_metrics = None
        self.target_metrics = None
        self.analysis_rows = []

        self._build_ui()

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        notebook = ttk.Notebook(self)
        notebook.grid(row=0, column=0, sticky="nsew", padx=8, pady=8)

        self.tab_input = ttk.Frame(notebook)
        self.tab_results = ttk.Frame(notebook)
        self.tab_help = ttk.Frame(notebook)

        notebook.add(self.tab_input, text="Input")
        notebook.add(self.tab_results, text="Analyse & Ergebnisse")
        notebook.add(self.tab_help, text="Hilfe")

        self._build_input_tab()
        self._build_results_tab()
        self._build_help_tab()

    def _build_input_tab(self):
        f = self.tab_input
        for i in range(2):
            f.columnconfigure(i, weight=1)
        for i in range(6):
            f.rowconfigure(i, weight=1)

        seq_frame = ttk.LabelFrame(f, text="Sequenz-Input")
        seq_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=6, pady=6)
        seq_frame.columnconfigure(0, weight=1)
        seq_frame.columnconfigure(1, weight=1)

        ttk.Label(seq_frame, text="Aminosäuresequenz").grid(row=0, column=0, sticky="w", padx=6, pady=(6, 2))
        ttk.Label(seq_frame, text="Coding-Sequenz im Source-Organismus").grid(row=0, column=1, sticky="w", padx=6, pady=(6, 2))

        self.aa_text = tk.Text(seq_frame, height=8, wrap="word")
        self.cds_text = tk.Text(seq_frame, height=8, wrap="word")
        self.aa_text.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        self.cds_text.grid(row=1, column=1, sticky="nsew", padx=6, pady=6)

        tables_frame = ttk.LabelFrame(f, text="Codonfrequenz-Tabellen")
        tables_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=6, pady=6)
        tables_frame.columnconfigure(0, weight=1)
        tables_frame.columnconfigure(1, weight=1)
        tables_frame.rowconfigure(1, weight=1)

        ttk.Label(tables_frame, text="Source-Organismus").grid(row=0, column=0, sticky="w", padx=6, pady=(6, 2))
        ttk.Label(tables_frame, text="Target-Organismus").grid(row=0, column=1, sticky="w", padx=6, pady=(6, 2))

        self.source_table_text = tk.Text(tables_frame, wrap="none")
        self.target_table_text = tk.Text(tables_frame, wrap="none")
        self.source_table_text.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        self.target_table_text.grid(row=1, column=1, sticky="nsew", padx=6, pady=6)

        controls = ttk.LabelFrame(f, text="Aktionen")
        controls.grid(row=2, column=0, columnspan=2, sticky="ew", padx=6, pady=6)
        for i in range(10):
            controls.columnconfigure(i, weight=1)

        ttk.Label(controls, text="Rare-Codon-Schwelle (RFN ≤)").grid(row=0, column=0, sticky="e", padx=6, pady=8)
        self.rare_var = tk.IntVar(value=1)
        ttk.Spinbox(controls, from_=1, to=6, textvariable=self.rare_var, width=5).grid(row=0, column=1, sticky="w", padx=6, pady=8)

        ttk.Button(controls, text="Beispieldaten laden", command=self.load_example).grid(row=0, column=2, padx=6, pady=8)
        ttk.Button(controls, text="Source-Tabelle aus Datei", command=lambda: self.load_table_into(self.source_table_text)).grid(row=0, column=3, padx=6, pady=8)
        ttk.Button(controls, text="Target-Tabelle aus Datei", command=lambda: self.load_table_into(self.target_table_text)).grid(row=0, column=4, padx=6, pady=8)
        ttk.Button(controls, text="Analyse starten", command=self.run_analysis).grid(row=0, column=5, padx=6, pady=8)
        ttk.Button(controls, text="Inputs leeren", command=self.clear_inputs).grid(row=0, column=6, padx=6, pady=8)

        status_frame = ttk.LabelFrame(f, text="Status")
        status_frame.grid(row=3, column=0, columnspan=2, sticky="ew", padx=6, pady=6)
        self.status_var = tk.StringVar(value="Bereit.")
        ttk.Label(status_frame, textvariable=self.status_var).grid(row=0, column=0, sticky="w", padx=8, pady=8)

    def _build_results_tab(self):
        f = self.tab_results
        f.columnconfigure(0, weight=1)
        f.rowconfigure(1, weight=1)
        f.rowconfigure(3, weight=1)

        summary = ttk.LabelFrame(f, text="Zusammenfassung")
        summary.grid(row=0, column=0, sticky="ew", padx=6, pady=6)
        summary.columnconfigure(0, weight=1)
        self.summary_text = tk.Text(summary, height=8, wrap="word")
        self.summary_text.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)

        table_frame = ttk.LabelFrame(f, text="Codon-für-Codon Analyse")
        table_frame.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        columns = [
            "Position", "AA", "Source Codon", "Source Frequency", "Source Score", "Source RFN",
            "Target Codon", "Target Frequency", "Target Score", "Target RFN", "Note"
        ]
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=20)
        for col in columns:
            self.tree.heading(col, text=col)
            width = 90
            if col in {"Note"}:
                width = 430
            elif col in {"Source Score", "Target Score"}:
                width = 100
            elif col in {"Source Frequency", "Target Frequency"}:
                width = 110
            self.tree.column(col, width=width, anchor="center")
        self.tree.grid(row=0, column=0, sticky="nsew")

        scroll_y = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        scroll_y.grid(row=0, column=1, sticky="ns")
        scroll_x = ttk.Scrollbar(table_frame, orient="horizontal", command=self.tree.xview)
        scroll_x.grid(row=1, column=0, sticky="ew")
        self.tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)

        export_frame = ttk.LabelFrame(f, text="Abgeleitete Zielsequenz")
        export_frame.grid(row=2, column=0, sticky="ew", padx=6, pady=6)
        export_frame.columnconfigure(0, weight=1)
        export_frame.columnconfigure(1, weight=1)

        ttk.Label(export_frame, text="Target Coding Sequence").grid(row=0, column=0, sticky="w", padx=6, pady=(6, 2))
        ttk.Label(export_frame, text="Target Amino Acid Sequence").grid(row=0, column=1, sticky="w", padx=6, pady=(6, 2))

        self.target_cds_text = tk.Text(export_frame, height=5, wrap="word")
        self.target_aa_text = tk.Text(export_frame, height=5, wrap="word")
        self.target_cds_text.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        self.target_aa_text.grid(row=1, column=1, sticky="nsew", padx=6, pady=6)

        btns = ttk.Frame(f)
        btns.grid(row=3, column=0, sticky="ew", padx=6, pady=6)
        ttk.Button(btns, text="Analyse als CSV exportieren", command=self.export_analysis_csv).pack(side="left", padx=4)
        ttk.Button(btns, text="Target-CDS speichern", command=self.export_target_cds).pack(side="left", padx=4)
        ttk.Button(btns, text="RFN-Tabellen speichern", command=self.export_rank_tables).pack(side="left", padx=4)

    def _build_help_tab(self):
        f = self.tab_help
        f.columnconfigure(0, weight=1)
        f.rowconfigure(0, weight=1)
        text = tk.Text(f, wrap="word")
        text.grid(row=0, column=0, sticky="nsew", padx=8, pady=8)
        text.insert("1.0", HELP_TEXT)
        text.configure(state="disabled")

    def set_status(self, text: str):
        self.status_var.set(text)
        self.update_idletasks()

    def clear_text(self, widget):
        widget.delete("1.0", "end")

    def load_example(self):
        self.clear_inputs()
        aa = "RKMF"
        cds = "CGTAAAATGTTT"
        self.aa_text.insert("1.0", aa)
        self.cds_text.insert("1.0", cds)
        self.source_table_text.insert("1.0", EXAMPLE_SOURCE)
        self.target_table_text.insert("1.0", EXAMPLE_TARGET)
        self.set_status("Beispieldaten geladen.")

    def clear_inputs(self):
        for widget in [self.aa_text, self.cds_text, self.source_table_text, self.target_table_text]:
            self.clear_text(widget)
        self.set_status("Inputs geleert.")

    def load_table_into(self, widget):
        path = filedialog.askopenfilename(
            title="Codonfrequenz-Tabelle auswählen",
            filetypes=[("Text/CSV/TSV", "*.txt *.csv *.tsv"), ("Alle Dateien", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                content = fh.read()
            self.clear_text(widget)
            widget.insert("1.0", content)
            self.set_status(f"Tabelle geladen: {path}")
        except Exception as exc:
            messagebox.showerror("Dateifehler", str(exc))

    def run_analysis(self):
        try:
            aa_input = clean_aa(self.aa_text.get("1.0", "end"))
            source_cds = clean_dna(self.cds_text.get("1.0", "end"))
            source_table = parse_frequency_table(self.source_table_text.get("1.0", "end"))
            target_table = parse_frequency_table(self.target_table_text.get("1.0", "end"))

            if not source_cds:
                raise ValueError("Bitte eine Coding-Sequenz eingeben.")
            if len(source_cds) % 3 != 0:
                raise ValueError("Die Coding-Sequenz muss durch 3 teilbar sein.")

            translated = translate_dna(source_cds)
            codon_count = len(source_cds) // 3

            if aa_input:
                if len(aa_input) != codon_count:
                    raise ValueError(
                        "Länge der Aminosäuresequenz passt nicht zur Coding-Sequenz.\n"
                        f"AA-Länge: {len(aa_input)}\n"
                        f"Codons in CDS: {codon_count}"
                    )
                desired_aa_seq = aa_input
            else:
                desired_aa_seq = translated

            self.source_metrics = compute_metrics(source_table)
            self.target_metrics = compute_metrics(target_table)

            rows, target_cds, target_aa, translated = analyze_sequence(
                source_cds=source_cds,
                desired_aa_seq=desired_aa_seq,
                source_metrics=self.source_metrics,
                target_metrics=self.target_metrics,
            )
            self.analysis_rows = rows

            self.populate_results(rows, target_cds, target_aa, translated, desired_aa_seq)
            self.set_status(f"Analyse abgeschlossen. {len(rows)} Codons verarbeitet.")
        except Exception as exc:
            messagebox.showerror("Analysefehler", str(exc))
            self.set_status("Analyse fehlgeschlagen.")

    def populate_results(self, rows, target_cds, target_aa, translated, desired_aa_seq):
        self.summary_text.delete("1.0", "end")
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.target_cds_text.delete("1.0", "end")
        self.target_aa_text.delete("1.0", "end")

        rare_count = sum(1 for r in rows if r["Source RFN"] <= self.rare_var.get())
        mismatch_count = sum(1 for i, aa in enumerate(desired_aa_seq) if translated[i] != aa)

        summary = (
            f"Anzahl Codons: {len(rows)}\n"
            f"Rare-Codons (RFN ≤ {self.rare_var.get()}): {rare_count}\n"
            f"AA-Mismatches zwischen Input-AA und Source-CDS: {mismatch_count}\n"
            f"Source-AA-Sequenz (aus CDS): {translated}\n"
            f"Soll-AA-Sequenz: {desired_aa_seq}\n"
            f"Abgeleitete Target-CDS: {target_cds}\n"
            f"Target-AA-Sequenz: {target_aa}\n\n"
            f"Definition:\n"
            f"- Score = Codonfrequenz / Summe der synonymen Codonfrequenzen\n"
            f"- RFN = Rang innerhalb einer Aminosäure, nach Score aufsteigend\n"
            f"- Zielmapping = Target-Codon mit möglichst ähnlichem Target-Score zum Source-Score,\n"
            f"  unter Beibehaltung der gewünschten Aminosäure\n"
        )
        self.summary_text.insert("1.0", summary)

        for row in rows:
            values = [
                row["Position"],
                row["AA"],
                row["Source Codon"],
                format_float(row["Source Frequency"]),
                format_float(row["Source Score"]),
                row["Source RFN"],
                row["Target Codon"],
                format_float(row["Target Frequency"]),
                format_float(row["Target Score"]),
                row["Target RFN"],
                row["Note"],
            ]
            tag = "rare" if row["Source RFN"] <= self.rare_var.get() else ""
            self.tree.insert("", "end", values=values, tags=(tag,))
        self.tree.tag_configure("rare", background="#ffe8e8")

        self.target_cds_text.insert("1.0", target_cds)
        self.target_aa_text.insert("1.0", target_aa)

    def export_analysis_csv(self):
        if not self.analysis_rows:
            messagebox.showinfo("Hinweis", "Noch keine Analyse vorhanden.")
            return
        path = filedialog.asksaveasfilename(
            title="Analyse speichern",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(self.analysis_rows[0].keys()))
            writer.writeheader()
            writer.writerows(self.analysis_rows)
        self.set_status(f"Analyse exportiert: {path}")

    def export_target_cds(self):
        seq = self.target_cds_text.get("1.0", "end").strip()
        if not seq:
            messagebox.showinfo("Hinweis", "Noch keine Target-CDS vorhanden.")
            return
        path = filedialog.asksaveasfilename(
            title="Target-CDS speichern",
            defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("FASTA", "*.fasta"), ("Alle Dateien", "*.*")],
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            if path.lower().endswith((".fasta", ".fa")):
                fh.write(">target_codon_sequence\n")
                fh.write(seq + "\n")
            else:
                fh.write(seq + "\n")
        self.set_status(f"Target-CDS gespeichert: {path}")

    def export_rank_tables(self):
        if self.source_metrics is None or self.target_metrics is None:
            messagebox.showinfo("Hinweis", "Noch keine RFN-Tabellen berechnet.")
            return
        path = filedialog.asksaveasfilename(
            title="RFN-Tabellen speichern",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not path:
            return
        rows = []
        for organism, metrics in [("source", self.source_metrics), ("target", self.target_metrics)]:
            for codon in sorted(metrics.keys()):
                info = metrics[codon]
                rows.append(
                    {
                        "organism": organism,
                        "codon": codon,
                        "aa": info["aa"],
                        "frequency": info["frequency"],
                        "score": info["score"],
                        "rfn": info["rfn"],
                        "synonymous_total": info["synonymous_total"],
                    }
                )
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        self.set_status(f"RFN-Tabellen gespeichert: {path}")

if __name__ == "__main__":
    app = RareCodonApp()
    app.mainloop()
