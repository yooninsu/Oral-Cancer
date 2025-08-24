import re
import csv
import argparse
import sys
import subprocess
from pathlib import Path

def read_s3_output_from_file(file_path):
    """파일에서 S3 출력 내용을 읽어옵니다."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except FileNotFoundError:
        print(f"파일을 찾을 수 없습니다: {file_path}")
        return None
    except Exception as e:
        print(f"파일 읽기 중 오류 발생: {e}")
        return None

def run_s3_command(bucket_path):
    """직접 AWS CLI 명령을 실행하여 S3 출력을 가져옵니다."""
    try:
        cmd = ['aws', 's3', 'ls', '--recursive', bucket_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"AWS CLI 명령 실행 실패: {e}")
        print(f"stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print("AWS CLI가 설치되지 않았거나 PATH에 없습니다.")
        return None

def read_from_stdin():
    """표준 입력에서 S3 출력을 읽어옵니다."""
    print("S3 출력을 붙여넣고 Ctrl+D (Linux/Mac) 또는 Ctrl+Z (Windows)를 눌러 입력을 완료하세요:")
    try:
        return sys.stdin.read()
    except KeyboardInterrupt:
        print("\n입력이 취소되었습니다.")
        return None

def create_manifest(s3_output, bucket_name="oral-cancer-ncc", output_csv_path='manifest.csv'):
    """
    S3 'ls' 명령어 출력을 파싱하여 manifest CSV 파일을 생성합니다.

    Args:
        s3_output (str): 'aws s3 ls --recursive' 명령어의 전체 출력 문자열.
        bucket_name (str): S3 버킷 이름.
        output_csv_path (str): 생성할 CSV 파일의 경로.
    """
    if not s3_output or not s3_output.strip():
        print("S3 출력 데이터가 비어있습니다.")
        return False

    # CSV 파일에 기록할 데이터를 저장할 리스트
    manifest_data = []
    
    # 샘플 ID를 추출하기 위한 정규표현식 (H로 시작하고 숫자 7자리)
    sample_id_pattern = re.compile(r'(H\d{7})')

    # S3 출력 문자열을 한 줄씩 처리
    for line_num, line in enumerate(s3_output.strip().split('\n'), 1):
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        
        # 파일 정보가 없는 줄은 건너뜁니다.
        if len(parts) < 4:
            continue
            
        file_path = ' '.join(parts[3:])  # 파일명에 공백이 있을 수 있음
        
        # 크기가 0이고 '/'로 끝나는 디렉터리 표시는 건너뜁니다.
        try:
            size = int(parts[2])
            if size == 0 and file_path.endswith('/'):
                continue
        except ValueError:
            print(f"Warning: 라인 {line_num}에서 파일 크기를 파싱할 수 없습니다: {line}")
            continue

        # 파일 경로에서 정보 추출
        path_parts = file_path.split('/')
        if len(path_parts) < 2:
            continue

        # 기본 정보 설정
        sample_id = 'N/A'
        data_type = 'N/A'
        group = 'N/A'
        read_pair = 'N/A'
        file_name = path_parts[-1]
        full_s3_path = f"s3://{bucket_name}/{file_path}"

        # 데이터 종류(data_type) 식별
        if len(path_parts) > 1:
            data_type = path_parts[1] # '16s', 'Olink', 'shotgun'
        
        # 샘플 ID 추출
        match = sample_id_pattern.search(file_name)
        if match:
            sample_id = match.group(1)

        # 데이터 종류에 따른 추가 정보 처리
        if data_type == '16s':
            if len(path_parts) > 2:
                group = path_parts[2] # 'Case' or 'Control'
        
        elif data_type == 'shotgun':
            if '_1.fastq.gz' in file_name:
                read_pair = 'R1'
            elif '_2.fastq.gz' in file_name:
                read_pair = 'R2'
        
        elif data_type == 'Olink':
            # Olink 데이터는 파일명에 샘플 ID가 없어 파일명을 기반으로 ID를 임시 지정합니다.
            if not match:
                sample_id = file_name.split('.')[0]

        # 추출된 정보를 딕셔너리로 저장
        row = {
            'sample_id': sample_id,
            'data_type': data_type,
            'group': group,
            'read_pair': read_pair,
            'file_name': file_name,
            's3_path': full_s3_path,
            'file_size': size
        }
        manifest_data.append(row)

    # CSV 파일 작성
    if not manifest_data:
        print("처리할 파일이 없습니다.")
        return False

    # CSV 헤더(컬럼명) 설정
    headers = ['sample_id', 'data_type', 'group', 'read_pair', 'file_name', 's3_path', 'file_size']
    
    try:
        with open(output_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(manifest_data)
            
        print(f"'{output_csv_path}' 파일이 성공적으로 생성되었습니다.")
        print(f"총 {len(manifest_data)}개의 파일 정보가 기록되었습니다.")
        return True
    except Exception as e:
        print(f"CSV 파일 작성 중 오류 발생: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='AWS S3 출력을 파싱하여 manifest CSV 파일을 생성합니다.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
사용 예시:
  # 파일에서 읽기
  python script.py -f s3_output.txt

  # AWS CLI로 직접 실행
  python script.py -s s3://your-bucket-name/

  # 표준 입력에서 읽기 (기본값)
  python script.py
  aws s3 ls --recursive s3://bucket-name/ | python script.py

  # 출력 파일명과 버킷명 지정
  python script.py -f input.txt -o output.csv -b my-bucket
        '''
    )
    
    parser.add_argument('-f', '--file', 
                        help='S3 출력이 저장된 파일 경로')
    parser.add_argument('-s', '--s3-path', 
                        help='직접 AWS CLI 명령을 실행할 S3 경로 (예: s3://bucket-name/)')
    parser.add_argument('-o', '--output', 
                        default='manifest.csv',
                        help='출력할 CSV 파일명 (기본값: manifest.csv)')
    parser.add_argument('-b', '--bucket', 
                        default='oral-cancer-ncc',
                        help='S3 버킷 이름 (기본값: oral-cancer-ncc)')
    
    args = parser.parse_args()
    
    # S3 출력 데이터 가져오기
    s3_output = None
    
    if args.file:
        print(f"파일에서 S3 출력 읽는 중: {args.file}")
        s3_output = read_s3_output_from_file(args.file)
    elif args.s3_path:
        print(f"AWS CLI 실행 중: aws s3 ls --recursive {args.s3_path}")
        s3_output = run_s3_command(args.s3_path)
    else:
        # 파이프라인 입력 확인
        if not sys.stdin.isatty():
            # 파이프된 입력이 있음
            s3_output = sys.stdin.read()
        else:
            # 대화형 입력
            s3_output = read_from_stdin()
    
    if s3_output is None:
        sys.exit(1)
    
    # Manifest 생성
    success = create_manifest(s3_output, args.bucket, args.output)
    
    if not success:
        sys.exit(1)
    
    print("\n처리 완료!")

if __name__ == "__main__":
    main()