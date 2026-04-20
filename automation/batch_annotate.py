#!/usr/bin/env python3
import argparse, sys, subprocess
from pathlib import Path
import pandas as pd
from datetime import datetime

PIPELINE_STEPS = [
    ("panelapp", "panelapp"), ("revel", "revel"), ("alphamissense", "alphamissense"),
    ("clinvar", "clinvar"), ("gnomad", "gnomad"), ("vep", "vep"), ("spliceai", "spliceai"),
]

STEP_TO_SCRIPT = {
    "panelapp": "annotate_panelapp_step.py", "revel": "annotate_revel_step.py",
    "alphamissense": "annotate_alphamissense_step.py", "clinvar": "annotate_clinvar_step.py",
    "gnomad": "annotate_gnomad_step.py", "vep": "create_small_vcf_for_vep.py",
    "spliceai": "annotate_spliceai_step.py",
}

STEP_TO_COLUMN = {
    "panelapp": "panelapp_result", "revel": "revel_result", "alphamissense": "alphamissense_result",
    "clinvar": "clinvar_result", "gnomad": "gnomad_result", "vep": "vep_input", "spliceai": "spliceai_result",
}

class BatchAnnotator:
    def __init__(self, sample_list_path, project_root=None):
        self.sample_list_path = Path(sample_list_path)
        if not self.sample_list_path.exists():
            raise FileNotFoundError(f"Not found: {sample_list_path}")
        self.project_root = self.sample_list_path.parent if project_root is None else Path(project_root)
        self.automation_dir = self.project_root / "automation"
        self.result_dir = self.project_root / "result"
        self.result_dir.mkdir(exist_ok=True)
        self.df = pd.read_csv(self.sample_list_path, sep="\t")
    
    def get_samples(self, sample_id=None, status=None):
        df = self.df.copy()
        if sample_id: df = df[df["sample_id"] == sample_id]
        if status: df = df[df["status"] == status]
        return df
    
    def update_status(self, sample_id, status, step=None):
        idx = self.df[self.df["sample_id"] == sample_id].index[0]
        self.df.at[idx, "status"] = status
        if step and step in STEP_TO_COLUMN:
            col = STEP_TO_COLUMN[step]
            if col in self.df.columns:
                self.df.at[idx, col] = "✓" if status == "completed" else ("✗" if status == "failed" else "→")
        self.df.to_csv(self.sample_list_path, sep="\t", index=False)
    
    def run_step(self, sample_id, step, step_input, step_output=None):
        if step == "vep": return True, "Manual VEP"
        script = self.automation_dir / STEP_TO_SCRIPT[step]
        if not script.exists(): return False, f"No script"
        try:
            cmd = [sys.executable, str(script), str(step_input)]
            if step_output: cmd.extend(["-o", str(step_output)])
            result = subprocess.run(cmd, cwd=str(self.project_root), capture_output=True, timeout=3600)
            return result.returncode == 0, "OK" if result.returncode == 0 else "Failed"
        except Exception as e: return False, str(e)
    
    def annotate_sample(self, sample_id, start_from="panelapp", skip_vep=True):
        sample = self.get_samples(sample_id=sample_id)
        if sample.empty: return False
        vcf_file = self.project_root / sample.iloc[0]["vcf_file"]
        if not vcf_file.exists(): return False
        self.update_status(sample_id, "in_progress")
        print(f"Processing: {sample_id}")
        start_idx = next((i for i, (s, _) in enumerate(PIPELINE_STEPS) if s == start_from), 0)
        current_input = vcf_file
        for step, _ in PIPELINE_STEPS[start_idx:]:
            if skip_vep and step == "vep": continue
            output_file = self.result_dir / f"{sample_id}_{step}.tsv"
            success, msg = self.run_step(sample_id, step, current_input, output_file)
            if not success: 
                self.update_status(sample_id, "failed", step)
                return False
            self.update_status(sample_id, "in_progress", step)
            if output_file.exists(): current_input = output_file
        self.update_status(sample_id, "completed")
        print(f"Done: {sample_id}")
        return True

def main():
    parser = argparse.ArgumentParser(description="Batch annotation")
    parser.add_argument("sample_list")
    parser.add_argument("--sample-id")
    parser.add_argument("--status", default="pending")
    parser.add_argument("--start-from", default="panelapp")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--reset", action="store_true")
    args = parser.parse_args()
    
    annotator = BatchAnnotator(args.sample_list)
    if args.list:
        print(annotator.df.to_string(index=False))
    elif args.reset:
        annotator.df["status"] = "pending"
        for col in STEP_TO_COLUMN.values():
            if col in annotator.df.columns: annotator.df[col] = ""
        annotator.df.to_csv(args.sample_list, sep="\t", index=False)
        print("Reset OK")
    elif args.sample_id:
        annotator.annotate_sample(args.sample_id, start_from=args.start_from)
    else:
        samples = annotator.get_samples(status=args.status)
        for _, sample in samples.iterrows():
            annotator.annotate_sample(sample["sample_id"], start_from=args.start_from)

if __name__ == "__main__":
    main()
