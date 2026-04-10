# ============================================================
# Mapping Coverage 전처리 & 샘플 병합 (Google Colab용) v3
# ============================================================
# 처리 내용:
#   1. ref_genome(C열)이 '-'인 gap 행 단순 삭제, A/T/G/C 행만 유지
#   2. _v1/_v2/_v3/_v4 등 모든 replicate를 동일 샘플 기준으로 A,C,G,T,total 합산 병합
#   3. 병합된 샘플을 각각 별도 시트로 저장
# ============================================================

from google.colab import files
import pandas as pd
from io import BytesIO
import re

# ── 파일 업로드 ──────────────────────────────────────────────
print("📂 Excel 파일을 업로드해 주세요...")
uploaded = files.upload()
filename = list(uploaded.keys())[0]
print(f"✅ 업로드 완료: {filename}\n")

# ── 전체 시트 읽기 ────────────────────────────────────────────
all_sheets = pd.read_excel(BytesIO(uploaded[filename]), sheet_name=None, header=0)
print(f"📋 발견된 시트 ({len(all_sheets)}개): {list(all_sheets.keys())}\n")

COL_NAMES = ['sample_id', 'number', 'ref_genome', 'A', 'C', 'G', 'T', 'total']
BASES = {'A', 'T', 'G', 'C'}
NUM_COLS = ['A', 'C', 'G', 'T']

def load_sheet(df):
    """시트 DataFrame을 표준 컬럼으로 변환"""
    df = df.copy()
    df.columns = COL_NAMES
    df['number'] = pd.to_numeric(df['number'], errors='coerce')
    df = df.dropna(subset=['number'])
    df['number'] = df['number'].astype(int)
    return df.reset_index(drop=True)

def deduplicate_by_atgc(df):
    """gap(-) 행 제거 → A/T/G/C 행만 유지, number 순 정렬"""
    return df[df['ref_genome'].isin(BASES)].sort_values('number').reset_index(drop=True)

def merge_replicates(dfs):
    """
    여러 replicate DataFrame(리스트)을 number 기준으로 outer join 후
    A, C, G, T, total 합산.
    ref_genome, sample_id는 첫 번째 df 기준 사용.
    """
    base = dfs[0][['number', 'ref_genome', 'sample_id', 'A', 'C', 'G', 'T', 'total']].copy()

    for i, df in enumerate(dfs[1:], start=2):
        base = pd.merge(
            base,
            df[['number', 'A', 'C', 'G', 'T', 'total']],
            on='number', how='outer', suffixes=('', f'_v{i}')
        ).fillna(0)
        for col in NUM_COLS + ['total']:
            base[col] = base[col].astype(int) + base[f'{col}_v{i}'].astype(int)
            base = base.drop(columns=[f'{col}_v{i}'])

    result = base[['sample_id', 'number', 'ref_genome', 'A', 'C', 'G', 'T', 'total']]
    return result.sort_values('number').reset_index(drop=True)

def base_name(name):
    """시트명에서 _v1, _v2, _v3, ... 제거"""
    return re.sub(r'_v\d+$', '', name, flags=re.IGNORECASE)

# ── v1/v2/v3/... 그룹 자동 탐지 ─────────────────────────────
groups = {}
for name in all_sheets:
    groups.setdefault(base_name(name), []).append(name)

print(f"🔗 병합 대상 샘플 그룹 ({len(groups)}개):")
for bn, sheets in groups.items():
    print(f"  {bn:35s} ← {sheets}")
print()

# ── 처리 및 결과 저장 ─────────────────────────────────────────
output = BytesIO()
with pd.ExcelWriter(output, engine='openpyxl') as writer:
    for bn, sheets in groups.items():
        sorted_sheets = sorted(sheets)

        dfs = [deduplicate_by_atgc(load_sheet(all_sheets[s])) for s in sorted_sheets]
        merged = merge_replicates(dfs)

        sheet_label = bn[:31]
        merged.to_excel(writer, sheet_name=sheet_label, index=False)

        counts = ' + '.join(str(len(d)) for d in dfs)
        print(f"✅ [{sheet_label}]  {counts} → 병합 {len(merged):,}개 위치  (replicate {len(dfs)}개)")

# ── 다운로드 ─────────────────────────────────────────────────
output.seek(0)
out_filename = filename.replace('.xlsx', '_merged.xlsx')
with open(out_filename, 'wb') as f:
    f.write(output.read())

files.download(out_filename)
print(f"\n🎉 완료! '{out_filename}' 파일 다운로드를 시작합니다.")
