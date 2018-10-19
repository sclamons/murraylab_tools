# -*- coding: utf-8 -*-
import io
import sys

__all__ = ['mt_open']

# Encoding nonsense (for detecting UTF-8 vs. UTF-16)
from codecs import BOM_UTF8, BOM_UTF16, BOM_UTF16_BE, BOM_UTF16_LE, \
                   BOM_UTF32, BOM_UTF32_BE, BOM_UTF32_LE

BOMS = (
    (BOM_UTF8, "UTF-8"),
    (BOM_UTF32, "UTF-32"),
    (BOM_UTF32_BE, "UTF-32-BE"),
    (BOM_UTF32_LE, "UTF-32-LE"),
    (BOM_UTF16, "UTF-16"),
    (BOM_UTF16_BE, "UTF-16-BE"),
    (BOM_UTF16_LE, "UTF-16-LE"),
)

def check_bom(data):
    return [encoding for bom, encoding in BOMS if data.startswith(bom)]


def mt_open(filename, setting_code):
    return_file = io.open(filename, setting_code,
                          encoding = "latin-1")# encoding)
    return return_file