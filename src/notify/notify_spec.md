# Module Specification: notify (Sink 멀티탭) (L3)

**Package**: [src/notify/](file:///home/ben/projects/ai4ref/src/notify/) | **Spec Version**: 1.0 | **Code Lines Constraint**: 각 모듈 <= 300줄 | **Spec Lines Constraint**: <= 200줄

## Summary

게이트(006)가 선별한 랜드마크를 features.yml 로 켜진 모든 출력 통로(sink)로 fan-out 하는 출력
레이어. 모든 통로는 동일 `Sink` 인터페이스('같은 플러그')를 구현한다. 켜진 통로만 활성화되고,
prereq(키·의존성) 미충족이면 자동 skip(fail-soft). 한 통로의 전송 실패는 격리되어 다른 통로·
다음 항목을 막지 않는다. 새 통로 추가 = 등록 한 줄, 본체(검색·선별) 불변(헌법 VIII).

## API Specifications

### `Sink` (base.py) — 통로 인터페이스
```python
class Sink(ABC):
    name: str
    def available(self) -> tuple[bool, str]: ...   # (사용가능, 사유). 미충족 → 멀티탭 skip
    @abstractmethod
    def send(self, item: LandmarkItem): ...         # 전달 후 참조(또는 None) 반환
```
* 본체는 통로 종류를 몰라도 동일하게 다룬다 (FR-001).

### `LandmarkItem` (base.py) — 전달 단위
* 필드: `kq, pmid, title, journal="", year="", reason="", url="", abstract=""`.
* `__post_init__`: `url` 미지정 시 pmid 에서 PubMed 링크 자동 파생 (FR-011).
* `key -> "kq:pmid"`: 중복 차단 식별자(§004 seen-set 연계).
* `abstract` = 게이트 판정용 **내부** 필드 — 렌더링에 쓰지 않는다 (FR-006).

### `build_sinks(features: dict) -> list[Sink]` (registry.py)
* `features["notify"]` 의 통로별 `enabled` 토글로 켜진 통로만 인스턴스화 (FR-002).
* 각 통로 `available()` 점검 → 미충족이면 `log.warning` 후 skip(활성 제외) (FR-003, fail-soft ①).
* `SINK_CLASSES = {"stdout": …, "telegram": …, "zotero": …}` — **새 통로 = 한 줄 등록**, 본체 0줄 (FR-009, 헌법 VIII).

### `fan_out(item, sinks: list) -> dict[str, ref|None]` (registry.py)
* 활성 통로 전체로 `item` 분배, 통로별 결과 `{name: 참조|None}` 반환 (FR-004).
* 통로 `send` 예외 → `out[name]=None`·`log.error`, **다른 통로·다음 항목 계속** (FR-005, fail-soft ②).

## 통로 구현 (sinks)

### `StdoutSink` (stdout_sink.py)
* 의존성 0·항상 동작 (개발·테스트 기본 구멍, FR-006). `available()=True`.
* `render(item) -> str`: 제목·메타(저널·연도)·한 줄 근거·PubMed 링크. 외부 통로가 **재사용**(FR-007).

### `TelegramSink` (telegram_sink.py)
* stdlib `urllib` 만(무비용, FR-010). `.env`: `TELEGRAM_BOT_TOKEN`·`TELEGRAM_CHAT_ID`.
* `available()`: 토큰/chat_id 부재 시 `(False, 사유)` → skip. 본문 = `render()` 재사용.

### `ZoteroSink` (zotero_sink.py)
* `pyzotero` **지연 로드**(미설치 시 skip, FR-010). `.env`: `ZOTERO_API_KEY`·`ZOTERO_USER_ID`(write).
* 지정 컬렉션(기본 `ai4ref-landmarks`) 생성/탐색 → 서지 항목 생성, `reason`→Extra 메모 보강 (FR-008).
* 메타·아이템 변환은 `collector.zotero_add`(fetch_meta·to_zotero_item) 재사용.

## 파이프라인 경계 (run.py)
* `run.py` 가 `build_sinks`/`fan_out` 호출. **안전장치**: 비-demo 실수집 + 게이트 stub 백엔드는
  모든 후보가 랜드마크가 되어 스팸 → 외부 통로 차단, stdout 만(§006/§007 경계, run.py 책임).

## 등가 기준 / 테스트
* 등가 = 기존 `alert.run` fan-out 동작 보존. 중복 차단(§004)·게이트(§006)·linked-file(§011)은 본 피처 밖.
* 단위테스트(`tests/unit/test_notify.py`): render(내부필드 제외)·LandmarkItem 파생·build_sinks 토글/
  skip·fan_out 전송 격리를 외부 호출 모의로 결정론 검증(텔레그램·Zotero 네트워크 미호출).
