"""Cache-first loading: checked-in cache, then the user cache, then live sources."""

from __future__ import annotations

import contextlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from snp_brief import pipeline
from snp_brief.pipeline import load_variant

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"


class LoadVariantTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)

    def test_checked_in_snp_reads_the_repo_cache_without_network(self) -> None:
        with mock.patch("snp_brief.pipeline.fetch_dbsnp", side_effect=AssertionError("network")):
            item = load_variant("rs6025", repo_cache=CACHE, user_cache=self.tmp)
        self.assertEqual(item.cache_dir, CACHE)
        self.assertEqual({call.origin for call in item.calls.values()}, {"cache"})

    def test_unknown_snp_is_fetched_live_into_the_user_cache(self) -> None:
        calls = []

        def fake_gather(rsid, *, cache_dir, refresh, offline):
            calls.append((rsid, cache_dir, refresh, offline))
            return mock.Mock(cache_dir=cache_dir)

        with mock.patch.object(pipeline, "gather", side_effect=fake_gather):
            load_variant("rs7903146", repo_cache=CACHE, user_cache=self.tmp)
        self.assertEqual(calls, [("rs7903146", self.tmp, False, False)])

    def test_user_cache_copy_wins_over_the_checked_in_one(self) -> None:
        shutil.copytree(CACHE / "rs6025", self.tmp / "rs6025")
        item = load_variant("rs6025", repo_cache=CACHE, user_cache=self.tmp)
        self.assertEqual(item.cache_dir, self.tmp)

    def test_refresh_refetches_into_the_user_cache_only(self) -> None:
        calls = []

        def fake_gather(rsid, *, cache_dir, refresh, offline):
            calls.append((cache_dir, refresh, offline))
            return mock.Mock(cache_dir=cache_dir)

        with mock.patch.object(pipeline, "gather", side_effect=fake_gather):
            load_variant("rs6025", repo_cache=CACHE, user_cache=self.tmp, refresh=True)
        self.assertEqual(calls, [(self.tmp, True, False)])

    def test_offline_unknown_snp_is_blank_not_an_error(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            item = load_variant("rs999999999", repo_cache=CACHE, user_cache=self.tmp, offline=True)
        self.assertEqual({call.origin for call in item.calls.values()}, {"no cache"})


if __name__ == "__main__":
    unittest.main()
