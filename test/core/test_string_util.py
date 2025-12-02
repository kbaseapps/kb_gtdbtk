import re
from datetime import datetime, timedelta
import time
from kb_gtdbtk.core.string_util import now_ISOish

def test_now_ISOish_format():
    """
    simple test for formatting
    """
    timestamp = now_ISOish()
    assert isinstance(timestamp, str)
    assert len(timestamp) == 15
    assert re.match(r"^\d{8}_\d{6}$", timestamp)

def test_now_ISOish_accuracy():
    """
    Generated timestamp should be accurate and match the stated format.
    """
    before_time = datetime.now()
    timestamp = now_ISOish()
    after_time = datetime.now()

    timestamp_time = datetime.strptime(timestamp, "%Y%m%d_%H%M%S")
    # Allow a small tolerance in case of second rollover
    assert before_time - timedelta(seconds=1) <= timestamp_time <= after_time + timedelta(seconds=1)

def test_now_ISOish_unique():
    """
    Generated timestamp should differ after a second.
    """
    first_time = now_ISOish()
    time.sleep(1)
    second_time = now_ISOish()
    assert first_time != second_time
