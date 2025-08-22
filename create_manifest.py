import pandas as pd
import os

# 1. AWS CLI로 생성한 텍스트 파일 이름
# (이전에 's3_all_samples_list.txt'로 저장했습니다)
s3_list_filename = 's3_all_samples_list.txt'

# 2. S3 버킷 이름 (본인 버킷 이름으로 수정)
bucket_name = 'labmember' 

manifest_data = []
try:
    with open(s3_list_filename, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"오류: '{s3_list_filename}' 파일을 찾을 수 없습니다.")
    print("먼저 aws s3 ls 명령어를 실행했는지 확인하세요.")
    exit()

for line in lines:
    parts = line.split()
    # 파일 경로가 있는 라인만 처리
    if len(parts) < 4 or not parts[3].endswith('.fastq.gz'):
        continue
    
    full_path = parts[3]
    filename = os.path.basename(full_path)
    
    # --- 핵심 로직: 파일명에서 샘플 ID와 그룹 정보 추출 ---
    
    # 1. 샘플 ID 추출: 파일명을 '_' 기준으로 나누고 첫 번째 요소(H1234567)를 가져옴
    sample_id = filename.split('_')[0]
    
    # 2. 그룹(Case/Control) 정보 추출: 전체 경로에 'Case'가 포함되어 있는지 확인
    group = 'Case' if '/Case/' in full_path else 'Control'
    
    # 3. S3 전체 경로(URI) 생성
    s3_uri = f's3://{bucket_name}/{full_path}'
    
    # ----------------------------------------------------
    
    manifest_data.append({
        'sample_id': sample_id,
        'group': group,
        's3_path': s3_uri,
        'filename': filename
    })

# DataFrame으로 변환 후 CSV로 저장
df = pd.DataFrame(manifest_data)
df.to_csv('sample_manifest.csv', index=False)

print("총 {len(df)}개의 샘플 정보를 처리하여 'sample_manifest.csv' 파일을 생성했습니다.")
print("파일 내용 미리보기:")
print(df.head())
