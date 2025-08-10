# zotero 모듈 리팩토링 계획

## 현재 문제점
- 235줄의 큰 단일 클래스
- 여러 책임이 혼재 (아이템, 컬렉션, 첨부, 배치처리)
- 테스트 및 유지보수 어려움

## 리팩토링 계획

### 1. 아이템 관련 (item_operations.py)
```python
class ZoteroItemManager:
    def build_item(self, paper_data)
    def search_by_pmid(self, pmid)
    def create_item(self, paper_data)
    def get_item(self, item_key)
```

### 2. 컬렉션 관련 (collection_operations.py)
```python
class ZoteroCollectionManager:
    def get_or_create_collection(self, name)
    def add_item_to_collection(self, item_key, collection_key)
    def list_collections(self)
```

### 3. 첨부 관리 (attachment_manager.py)
```python
class ZoteroAttachmentManager:
    def attach_pdf(self, item_key, pdf_path, pmid)
    def check_pdf_exists(self, pdf_path)
```

### 4. 중복 감지 (duplicate_detector.py)
```python
class ZoteroDuplicateDetector:
    def detect_existing_item(self, pmid)
    def resolve_conflicts(self, existing_key, new_data)
```

### 5. 배치 처리 (batch_processor.py)
```python
class ZoteroBatchProcessor:
    def process_paper_collections(self, paper_id)
    def upload_all_papers(self, batch_size)
    def show_sync_status(self)
```

### 6. 통합 클라이언트 (zotero_client.py)
```python
class ZoteroClient:
    def __init__(self):
        self.item_manager = ZoteroItemManager()
        self.collection_manager = ZoteroCollectionManager()
        self.attachment_manager = ZoteroAttachmentManager()
        self.duplicate_detector = ZoteroDuplicateDetector()
        self.batch_processor = ZoteroBatchProcessor()
```

## 장점
1. 단일 책임 원칙 준수
2. 테스트 용이성 향상
3. 코드 재사용성 증대
4. 유지보수성 개선
5. 협업 시 충돌 감소

## 단점
1. 초기 설정 복잡도 증가
2. 임포트 관리 필요
3. 파일 수 증가

## 업계 트렌드 부합도
✅ 현대적 Python 프로젝트 구조
✅ 마이크로서비스 아키텍처 준비
✅ 테스트 주도 개발 지원
✅ 타입 힌트 적용 용이
