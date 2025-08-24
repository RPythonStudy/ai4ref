import psycopg2
from common.database import get_db_connection
from common.logger import log_debug, log_error

def show_table_structure():
    """
    생성된 테이블들의 구조를 표 형태로 보여줍니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()
        
        # 테이블 목록 가져오기 (실제 존재하는 스키마 테이블들)
        tables = ['papers', 'collections', 'collection_pmid', 'unique_pmid', 'filtered_pmid', 'paper_collection']
        
        for table_name in tables:
            # 테이블 구조 조회
            cur.execute("""
                SELECT 
                    column_name,
                    data_type,
                    character_maximum_length,
                    is_nullable,
                    column_default
                FROM information_schema.columns 
                WHERE table_name = %s
                ORDER BY ordinal_position;
            """, (table_name,))
            
            columns = cur.fetchall()
            
            if not columns:
                print(f"{table_name} 테이블이 존재하지 않습니다")
                continue
            
            # 테이블 헤더
            print("\n" + "="*80)
            print(f"{table_name.upper()} 테이블 구조".center(80))
            print("="*80)
            print(f"{'컬럼명':<20} {'타입':<15} {'길이':<8} {'NULL':<8} {'기본값':<20}")
            print("-"*80)
            
            # 테이블 데이터
            for col in columns:
                column_name = col[0]
                data_type = col[1]
                max_length = str(col[2]) if col[2] else ""
                is_nullable = "YES" if col[3] == "YES" else "NO"
                default_val = str(col[4])[:18] + "..." if col[4] and len(str(col[4])) > 18 else str(col[4]) if col[4] else ""
                
                print(f"{column_name:<20} {data_type:<15} {max_length:<8} {is_nullable:<8} {default_val:<20}")
            
            # 인덱스 정보
            cur.execute("""
                SELECT indexname, indexdef 
                FROM pg_indexes 
                WHERE tablename = %s
                ORDER BY indexname;
            """, (table_name,))
            
            indexes = cur.fetchall()
            
            if indexes:
                print("\n" + "-"*80)
                print("인덱스 정보".center(80))
                print("-"*80)
                for idx_name, idx_def in indexes:
                    print(f"- {idx_name}")
                    if "UNIQUE" in idx_def:
                        print(f"  (UNIQUE 제약조건)")
                    print()
            
            print("="*80 + "\n")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"테이블 구조 조회 중 오류: {e}")

def get_database_summary():
    """
    데이터베이스의 간단한 요약 정보를 반환합니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()
        
        summary = {}
        tables = ['papers', 'collections', 'collection_pmid', 'unique_pmid', 'filtered_pmid', 'paper_collection']
        
        for table in tables:
            cur.execute(f"SELECT COUNT(*) FROM {table};")
            count = cur.fetchone()[0]
            summary[table] = count
        
        cur.close()
        conn.close()
        return summary
        
    except psycopg2.Error as e:
        log_error(f"데이터베이스 요약 조회 중 오류: {e}")
        return {}

def show_table_columns():
    """
    각 테이블의 컬럼명만 간단히 보여줍니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()
        
        # 테이블 목록 가져오기 (실제 존재하는 스키마 테이블들)
        tables = ['papers', 'collections', 'collection_pmid', 'unique_pmid', 'filtered_pmid', 'paper_collection']
        
        print("\n" + "="*50)
        print("데이터베이스 테이블 구조 (컬럼명)".center(50))
        print("="*50)
        
        for table_name in tables:
            # 테이블 컬럼 조회
            cur.execute("""
                SELECT column_name, data_type
                FROM information_schema.columns 
                WHERE table_name = %s
                ORDER BY ordinal_position;
            """, (table_name,))
            
            columns = cur.fetchall()
            
            if not columns:
                print(f"\n{table_name.upper()}: 테이블이 존재하지 않습니다")
                continue
            
            print(f"\n📋 {table_name.upper()} ({len(columns)}개 컬럼)")
            print("-" * 30)
            
            # 컬럼명과 타입을 간단히 표시
            for col_name, col_type in columns:
                # 타입 간략화
                simple_type = col_type.split('(')[0]  # varchar(255) -> varchar
                if simple_type == 'character varying':
                    simple_type = 'varchar'
                elif simple_type == 'timestamp without time zone':
                    simple_type = 'timestamp'
                
                print(f"  • {col_name:<20} ({simple_type})")
        
        print("\n" + "="*50 + "\n")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"테이블 컬럼 조회 중 오류: {e}")

def show_database_summary():
    """
    데이터베이스 주요 지표와 컬럼별 데이터 수를 보기 좋게 출력합니다.
    """
    try:
        conn = psycopg2.connect(
            host="localhost",
            database="ai4ref",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()
        
        print("\n" + "="*60)
        print("데이터베이스 상태 요약".center(60))
        print("="*60)
        
        # 1. PAPERS 테이블 주요 지표
        cur.execute("SELECT COUNT(*) FROM papers;")
        total_papers = cur.fetchone()[0]
        
        # 데이터 유무와 관계없이 모든 컬럼 통계 조회
        cur.execute("SELECT COUNT(*) FROM papers WHERE doi IS NOT NULL AND TRIM(doi) != '';")
        papers_with_doi = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM papers WHERE pmid IS NOT NULL AND TRIM(pmid) != '';")
        papers_with_pmid = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM papers WHERE title IS NOT NULL AND TRIM(title) != '';")
        papers_with_title = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM papers WHERE pdf_url IS NOT NULL AND TRIM(pdf_url) != '';")
        papers_with_pdf = cur.fetchone()[0]
        
        print(f"📄 PAPERS: 총 {total_papers:,}개")
        if total_papers > 0:
            print(f"   ├─ DOI 보유: {papers_with_doi:,}개 ({papers_with_doi/total_papers*100:.1f}%)")
            print(f"   ├─ PMID 보유: {papers_with_pmid:,}개 ({papers_with_pmid/total_papers*100:.1f}%)")
            print(f"   ├─ 제목 보유: {papers_with_title:,}개 ({papers_with_title/total_papers*100:.1f}%)")
            print(f"   └─ PDF 링크: {papers_with_pdf:,}개 ({papers_with_pdf/total_papers*100:.1f}%)")
        else:
            print(f"   ├─ DOI 보유: {papers_with_doi:,}개 (0.0%)")
            print(f"   ├─ PMID 보유: {papers_with_pmid:,}개 (0.0%)")
            print(f"   ├─ 제목 보유: {papers_with_title:,}개 (0.0%)")
            print(f"   └─ PDF 링크: {papers_with_pdf:,}개 (0.0%)")
        
        # 2. COLLECTIONS 테이블 주요 지표
        cur.execute("SELECT COUNT(*) FROM collections;")
        total_collections = cur.fetchone()[0]
        
        # 데이터 유무와 관계없이 모든 컬럼 통계 조회
        cur.execute("SELECT COUNT(*) FROM collections WHERE enabled = true;")
        enabled_collections = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM collections WHERE parent_id IS NULL;")
        root_collections = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM collections WHERE zotero_key IS NOT NULL AND TRIM(zotero_key) != '';")
        zotero_collections = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM collections WHERE search_term IS NOT NULL AND TRIM(search_term) != '';")
        searchable_collections = cur.fetchone()[0]
        
        print(f"📁 COLLECTIONS: 총 {total_collections:,}개")
        if total_collections > 0:
            print(f"   ├─ 활성화됨: {enabled_collections:,}개 ({enabled_collections/total_collections*100:.1f}%)")
            print(f"   ├─ 최상위: {root_collections:,}개")
            print(f"   ├─ Zotero 연동: {zotero_collections:,}개 ({zotero_collections/total_collections*100:.1f}%)")
            print(f"   └─ 검색 설정: {searchable_collections:,}개 ({searchable_collections/total_collections*100:.1f}%)")
        else:
            print(f"   ├─ 활성화됨: {enabled_collections:,}개 (0.0%)")
            print(f"   ├─ 최상위: {root_collections:,}개")
            print(f"   ├─ Zotero 연동: {zotero_collections:,}개 (0.0%)")
            print(f"   └─ 검색 설정: {searchable_collections:,}개 (0.0%)")
        
        # 3. PAPER_COLLECTION 관계 테이블
        cur.execute("SELECT COUNT(*) FROM paper_collection;")
        total_relations = cur.fetchone()[0]
        
        # 데이터 유무와 관계없이 모든 통계 조회
        cur.execute("SELECT COUNT(DISTINCT paper_id) FROM paper_collection;")
        unique_papers = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT collection_id) FROM paper_collection;")
        unique_collections = cur.fetchone()[0]
        
        print(f"🔗 PAPER_COLLECTION: 총 {total_relations:,}개 관계")
        if total_relations > 0:
            print(f"   ├─ 연결된 논문: {unique_papers:,}개")
            print(f"   └─ 연결된 컬렉션: {unique_collections:,}개")
        else:
            print(f"   ├─ 연결된 논문: {unique_papers:,}개")
            print(f"   └─ 연결된 컬렉션: {unique_collections:,}개")
        
        # 4. 워크플로우 테이블들 상태
        cur.execute("SELECT COUNT(*) FROM collection_pmid;")
        collection_pmid_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM unique_pmid;")
        unique_pmid_count = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(*) FROM filtered_pmid;")
        filtered_pmid_count = cur.fetchone()[0]
        
        print(f"⚙️  워크플로우 상태")
        print(f"   ├─ collection_pmid: {collection_pmid_count:,}개 (Step 2: 검색 결과)")
        print(f"   ├─ unique_pmid: {unique_pmid_count:,}개 (Step 3a: 중복 제거)")
        print(f"   └─ filtered_pmid: {filtered_pmid_count:,}개 (Step 3b: 최종 필터링)")
        
        print("="*60 + "\n")
        
        cur.close()
        conn.close()
        
    except psycopg2.Error as e:
        log_error(f"데이터베이스 요약 조회 중 오류: {e}")
        print("\n데이터베이스 연결 오류로 요약 정보를 가져올 수 없습니다.\n")

def check_id_duplicates() -> None:
    """ID 중복 검사 - 반복구조"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_debug("=== ID 중복 검사 ===")
        
        # 검사할 테이블과 컬럼 설정
        checks = [
            {
                'table': 'papers',
                'column': 'pmid',
                'label': 'PMID',
                'condition': "pmid IS NOT NULL AND TRIM(pmid) != ''"
            },
            {
                'table': 'papers', 
                'column': 'doi',
                'label': 'DOI',
                'condition': "doi IS NOT NULL AND TRIM(doi) != ''"
            },
            {
                'table': 'papers', 
                'column': 'pmcid',
                'label': 'PMCID',
                'condition': "pmcid IS NOT NULL AND TRIM(pmcid) != ''"
            },
            {
                'table': 'papers', 
                'column': 'arxiv_id',
                'label': 'arXiv ID',
                'condition': "arxiv_id IS NOT NULL AND TRIM(arxiv_id) != ''"
            },
            {
                'table': 'paper_collection',
                'column': 'paper_id', 
                'label': 'Paper ID',
                'condition': "paper_id IS NOT NULL"
            },
            {
                'table': 'collections',
                'column': 'zotero_key',
                'label': 'Zotero Key',
                'condition': "zotero_key IS NOT NULL AND TRIM(zotero_key) != ''"
            },
            {
                'table': 'collection_pmid',
                'column': 'pmid',
                'label': 'Collection PMID',
                'condition': "pmid IS NOT NULL AND TRIM(pmid) != ''"
            },
            {
                'table': 'unique_pmid',
                'column': 'pmid',
                'label': 'Unique PMID',
                'condition': "pmid IS NOT NULL AND TRIM(pmid) != ''"
            },
            {
                'table': 'filtered_pmid',
                'column': 'pmid',
                'label': 'Filtered PMID',
                'condition': "pmid IS NOT NULL AND TRIM(pmid) != ''"
            }
        ]
        
        for i, check in enumerate(checks):
            try:
                if i > 0:
                    log_debug("")
                
                table_name = check['table'].upper()
                log_debug(f"[{table_name}] 테이블 - {check['label']} 중복 검사")
                
                cur.execute(f"""
                    SELECT {check['column']}, COUNT(*) as duplicate_count
                    FROM {check['table']}
                    WHERE {check['condition']}
                    GROUP BY {check['column']}
                    HAVING COUNT(*) > 1
                    ORDER BY COUNT(*) DESC, {check['column']}
                    LIMIT 10;
                """)
                
                duplicates = cur.fetchall()
                if duplicates:
                    log_debug(f"  {check['label']} 중복 발견: {len(duplicates)}개 (상위 10개 표시)")
                    for value, count in duplicates:
                        log_debug(f"    - {value}: {count}회 중복")
                else:
                    log_debug(f"  {check['label']} 중복 없음")
                    
            except Exception as e:
                log_debug(f"  {check['table']}.{check['column']} 검사 오류: {e}")
        
        log_debug("")
        
    finally:
        cur.close()
        conn.close()

def check_column_stats() -> None:
    """컬럼별 통계 출력"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_debug("=== 컬럼 통계 ===")
        
        # 1. 모든 테이블 목록 조회
        cur.execute("""
            SELECT tablename 
            FROM pg_tables 
            WHERE schemaname = 'public'
            ORDER BY tablename;
        """)
        tables = [row[0] for row in cur.fetchall()]
        
        for table_name in tables:
            # 2. 각 테이블의 컬럼 정보 조회
            cur.execute("""
                SELECT column_name, data_type
                FROM information_schema.columns 
                WHERE table_schema = 'public' 
                AND table_name = %s
                ORDER BY ordinal_position;
            """, (table_name,))
            columns = cur.fetchall()
            
            # 3. 테이블 총 행 수 조회
            cur.execute(f"SELECT COUNT(*) FROM {table_name};")
            result = cur.fetchone()
            total_rows = result[0] if result else 0
            
            log_debug(f"")
            log_debug(f"[{table_name.upper()}] 총 {total_rows:,}개 행")
            
            # 4. 각 컬럼별 데이터 개수 조회 및 출력
            for column_name, data_type in columns:
                if data_type in ['text', 'character varying', 'varchar']:
                    # 텍스트 컬럼: NULL이 아니고 빈 문자열이 아닌 경우
                    cur.execute(f"""
                        SELECT COUNT(*) 
                        FROM {table_name} 
                        WHERE {column_name} IS NOT NULL AND TRIM({column_name}) != '';
                    """)
                else:
                    # 숫자, 날짜 등: NULL이 아닌 경우
                    cur.execute(f"""
                        SELECT COUNT(*) 
                        FROM {table_name} 
                        WHERE {column_name} IS NOT NULL;
                    """)
                
                result = cur.fetchone()
                count = result[0] if result else 0
                percentage = (count / total_rows * 100) if total_rows > 0 else 0
                log_debug(f"  - {column_name} ({data_type}): {count:,}개 ({percentage:.1f}%)")
        
        log_debug("")
        
    finally:
        cur.close()
        conn.close()

def check_workflow_status() -> None:
    """새로운 논문 수집 워크플로우 상태 분석"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_debug("=== 워크플로우 상태 분석 ===")
        
        # Step 2: collection_pmid 상태 분석
        cur.execute("SELECT COUNT(*) FROM collection_pmid;")
        collection_pmid_total = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT collection_id) FROM collection_pmid;")
        collection_pmid_collections = cur.fetchone()[0]
        
        log_debug(f"[Step 2] collection_pmid: {collection_pmid_total:,}개 PMID")
        log_debug(f"  - 관련 컬렉션: {collection_pmid_collections:,}개")
        
        # Step 3a: unique_pmid 상태 분석
        cur.execute("SELECT COUNT(*) FROM unique_pmid;")
        unique_pmid_total = cur.fetchone()[0]
        
        cur.execute("SELECT COUNT(DISTINCT first_collection_id) FROM unique_pmid;")
        unique_pmid_collections = cur.fetchone()[0]
        
        log_debug(f"[Step 3a] unique_pmid: {unique_pmid_total:,}개 중복제거 PMID")
        log_debug(f"  - 관련 컬렉션: {unique_pmid_collections:,}개")
        
        # Step 3b: filtered_pmid 상태 분석
        cur.execute("SELECT COUNT(*) FROM filtered_pmid;")
        filtered_pmid_total = cur.fetchone()[0]
        
        log_debug(f"[Step 3b] filtered_pmid: {filtered_pmid_total:,}개 최종필터링 PMID")
        
        # 컬렉션별 워크플로우 진행 상태
        cur.execute("""
            SELECT 
                c.name as collection_name,
                COUNT(cp.pmid) as step2_count,
                COUNT(up.pmid) as step3a_count,
                COUNT(fp.pmid) as step3b_count
            FROM collections c
            LEFT JOIN collection_pmid cp ON c.id = cp.collection_id
            LEFT JOIN unique_pmid up ON c.id = up.first_collection_id
            LEFT JOIN filtered_pmid fp ON fp.pmid = up.pmid
            GROUP BY c.id, c.name
            HAVING COUNT(cp.pmid) > 0 OR COUNT(up.pmid) > 0 OR COUNT(fp.pmid) > 0
            ORDER BY step2_count DESC
            LIMIT 10;
        """)
        
        results = cur.fetchall()
        if results:
            log_debug("")
            log_debug("[컬렉션별] 워크플로우 진행 상태 (상위 10개)")
            for name, step2, step3a, step3b in results:
                log_debug(f"  - {name}: Step2={step2}, Step3a={step3a}, Step3b={step3b}")
        
        log_debug("")
        
    except Exception as e:
        log_debug(f"워크플로우 상태 분석 오류: {e}")
    finally:
        cur.close()
        conn.close()

def check_collection_hierarchy() -> None:
    """컬렉션 계층 구조 분석"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_debug("=== 컬렉션 계층 구조 ===")
        
        # 최상위 컬렉션 (parent_id가 NULL)
        cur.execute("""
            SELECT COUNT(*) 
            FROM collections 
            WHERE parent_id IS NULL;
        """)
        root_count = cur.fetchone()[0]
        log_debug(f"최상위 컬렉션: {root_count}개")
        
        # 하위 컬렉션이 있는 컬렉션
        cur.execute("""
            SELECT COUNT(DISTINCT parent_id) 
            FROM collections 
            WHERE parent_id IS NOT NULL;
        """)
        parent_count = cur.fetchone()[0]
        log_debug(f"하위 컬렉션을 가진 컬렉션: {parent_count}개")
        
        # 최대 깊이 계산
        cur.execute("""
            WITH RECURSIVE collection_tree AS (
                SELECT id, name, parent_id, 1 as level
                FROM collections 
                WHERE parent_id IS NULL
                
                UNION ALL
                
                SELECT c.id, c.name, c.parent_id, ct.level + 1
                FROM collections c
                JOIN collection_tree ct ON c.parent_id = ct.id
            )
            SELECT MAX(level) as max_depth
            FROM collection_tree;
        """)
        
        result = cur.fetchone()
        max_depth = result[0] if result and result[0] else 0
        log_debug(f"최대 계층 깊이: {max_depth}단계")
        
        # 계층별 분포
        cur.execute("""
            WITH RECURSIVE collection_tree AS (
                SELECT id, name, parent_id, 1 as level
                FROM collections 
                WHERE parent_id IS NULL
                
                UNION ALL
                
                SELECT c.id, c.name, c.parent_id, ct.level + 1
                FROM collections c
                JOIN collection_tree ct ON c.parent_id = ct.id
            )
            SELECT level, COUNT(*) as count
            FROM collection_tree
            GROUP BY level
            ORDER BY level;
        """)
        
        results = cur.fetchall()
        if results:
            log_debug("계층별 분포:")
            for level, count in results:
                log_debug(f"  - {level}단계: {count}개")
        
        log_debug("")
        
    except Exception as e:
        log_debug(f"컬렉션 계층 구조 분석 오류: {e}")
    finally:
        cur.close()
        conn.close()

def check_cross_table_consistency() -> None:
    """테이블 간 데이터 일관성 검사"""
    conn = get_db_connection()
    cur = conn.cursor()
    
    try:
        log_debug("=== 테이블 간 일관성 검사 ===")
        
        # paper_collection의 paper_id가 papers에 존재하는지 확인
        cur.execute("""
            SELECT COUNT(*) 
            FROM paper_collection pc
            LEFT JOIN papers p ON pc.paper_id = p.id
            WHERE p.id IS NULL;
        """)
        orphan_papers = cur.fetchone()[0]
        log_debug(f"[paper_collection] 존재하지 않는 paper_id 참조: {orphan_papers}개")
        
        # paper_collection의 collection_id가 collections에 존재하는지 확인
        cur.execute("""
            SELECT COUNT(*) 
            FROM paper_collection pc
            LEFT JOIN collections c ON pc.collection_id = c.id
            WHERE c.id IS NULL;
        """)
        orphan_collections = cur.fetchone()[0]
        log_debug(f"[paper_collection] 존재하지 않는 collection_id 참조: {orphan_collections}개")
        
        # collection_pmid의 collection_id가 collections에 존재하는지 확인
        cur.execute("""
            SELECT COUNT(*) 
            FROM collection_pmid cp
            LEFT JOIN collections c ON cp.collection_id = c.id
            WHERE c.id IS NULL;
        """)
        orphan_collection_pmid = cur.fetchone()[0]
        log_debug(f"[collection_pmid] 존재하지 않는 collection_id 참조: {orphan_collection_pmid}개")
        
        # unique_pmid의 first_collection_id가 collections에 존재하는지 확인
        cur.execute("""
            SELECT COUNT(*) 
            FROM unique_pmid up
            LEFT JOIN collections c ON up.first_collection_id = c.id
            WHERE c.id IS NULL;
        """)
        orphan_unique_pmid = cur.fetchone()[0]
        log_debug(f"[unique_pmid] 존재하지 않는 first_collection_id 참조: {orphan_unique_pmid}개")
        
        # collections의 parent_id 순환 참조 검사
        cur.execute("""
            WITH RECURSIVE cycle_check AS (
                SELECT id, parent_id, ARRAY[id] as path
                FROM collections
                WHERE parent_id IS NOT NULL
                
                UNION ALL
                
                SELECT c.id, c.parent_id, cc.path || c.id
                FROM collections c
                JOIN cycle_check cc ON c.parent_id = cc.id
                WHERE c.id != ALL(cc.path)
            )
            SELECT COUNT(*)
            FROM cycle_check cc
            JOIN collections c ON cc.id = c.id
            WHERE c.id = ANY(cc.path[1:array_length(cc.path, 1)-1]);
        """)
        
        cycle_count = cur.fetchone()[0]
        log_debug(f"[collections] 순환 참조 발견: {cycle_count}개")
        
        # 중복 조합 검사
        cur.execute("""
            SELECT COUNT(*) 
            FROM (
                SELECT paper_id, collection_id, COUNT(*)
                FROM paper_collection
                GROUP BY paper_id, collection_id
                HAVING COUNT(*) > 1
            ) duplicates;
        """)
        
        duplicate_relations = cur.fetchone()[0]
        log_debug(f"[paper_collection] 중복 관계: {duplicate_relations}개")
        
        # collection_pmid 중복 검사
        cur.execute("""
            SELECT COUNT(*) 
            FROM (
                SELECT collection_id, pmid, COUNT(*)
                FROM collection_pmid
                GROUP BY collection_id, pmid
                HAVING COUNT(*) > 1
            ) duplicates;
        """)
        
        duplicate_collection_pmid = cur.fetchone()[0]
        log_debug(f"[collection_pmid] 중복 항목: {duplicate_collection_pmid}개")
        
        # unique_pmid 중복 검사 - pmid 컬럼이 PRIMARY KEY이므로 first_collection_id 기준으로 검사
        cur.execute("""
            SELECT COUNT(*) 
            FROM (
                SELECT first_collection_id, pmid, COUNT(*)
                FROM unique_pmid
                GROUP BY first_collection_id, pmid
                HAVING COUNT(*) > 1
            ) duplicates;
        """)
        
        duplicate_unique_pmid = cur.fetchone()[0]
        log_debug(f"[unique_pmid] 중복 항목: {duplicate_unique_pmid}개")
        
        # filtered_pmid 중복 검사
        cur.execute("""
            SELECT COUNT(*) 
            FROM (
                SELECT pmid, COUNT(*)
                FROM filtered_pmid
                GROUP BY pmid
                HAVING COUNT(*) > 1
            ) duplicates;
        """)
        
        duplicate_filtered_pmid = cur.fetchone()[0]
        log_debug(f"[filtered_pmid] 중복 항목: {duplicate_filtered_pmid}개")
        
        log_debug("")
        
    except Exception as e:
        log_debug(f"테이블 간 일관성 검사 오류: {e}")
    finally:
        cur.close()
        conn.close()

def main():
    """
    메인 함수 - 사용자가 원하는 분석을 선택할 수 있습니다.
    """
    import sys
    
    if len(sys.argv) < 2:
        print("\n사용법:")
        print("  python analyze_database.py summary     # 데이터베이스 상태 요약")
        print("  python analyze_database.py columns     # 테이블 컬럼 구조")
        print("  python analyze_database.py structure   # 상세 테이블 구조")
        print("  python analyze_database.py stats       # 컬럼별 통계")
        print("  python analyze_database.py duplicates  # ID 중복 검사")
        print("  python analyze_database.py workflow    # 워크플로우 상태")
        print("  python analyze_database.py hierarchy   # 컬렉션 계층")
        print("  python analyze_database.py consistency # 테이블 간 일관성")
        print("  python analyze_database.py all         # 전체 분석")
        return
    
    command = sys.argv[1].lower()
    
    if command == "summary":
        show_database_summary()
    elif command == "columns":
        show_table_columns()
    elif command == "structure":
        show_table_structure()
    elif command == "stats":
        check_column_stats()
    elif command == "duplicates":
        check_id_duplicates()
    elif command == "workflow":
        check_workflow_status()
    elif command == "hierarchy":
        check_collection_hierarchy()
    elif command == "consistency":
        check_cross_table_consistency()
    elif command == "all":
        print("=== 전체 데이터베이스 분석 시작 ===\n")
        show_database_summary()
        show_table_columns()
        check_column_stats()
        check_id_duplicates()
        check_workflow_status()
        check_collection_hierarchy()
        check_cross_table_consistency()
        print("=== 전체 분석 완료 ===")
    else:
        print(f"알 수 없는 명령어: {command}")
        print("사용 가능한 명령어: summary, columns, structure, stats, duplicates, workflow, hierarchy, consistency, all")

if __name__ == "__main__":
    main()
