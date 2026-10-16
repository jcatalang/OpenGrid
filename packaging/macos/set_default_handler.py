"""Make a bundle the default application of a content type (Launch Services), so the Finder
opens those files with it even when other bundles claim the type:
python set_default_handler.py org.opengrid.project org.opengrid.OpenGrid.dev"""

import ctypes
import ctypes.util
import sys

K_LS_ROLES_ALL = 0xFFFFFFFF
UTF8 = 0x08000100  # kCFStringEncodingUTF8


def main() -> int:
    content_type, bundle_id = sys.argv[1], sys.argv[2]
    core = ctypes.cdll.LoadLibrary(ctypes.util.find_library("CoreFoundation"))
    services = ctypes.cdll.LoadLibrary(ctypes.util.find_library("CoreServices"))
    core.CFStringCreateWithCString.restype = ctypes.c_void_p
    core.CFStringCreateWithCString.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint32]
    core.CFRelease.argtypes = [ctypes.c_void_p]
    services.LSSetDefaultRoleHandlerForContentType.restype = ctypes.c_int32
    services.LSSetDefaultRoleHandlerForContentType.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_void_p,
    ]
    uti = core.CFStringCreateWithCString(None, content_type.encode(), UTF8)
    bundle = core.CFStringCreateWithCString(None, bundle_id.encode(), UTF8)
    status = services.LSSetDefaultRoleHandlerForContentType(uti, K_LS_ROLES_ALL, bundle)
    core.CFRelease(uti)
    core.CFRelease(bundle)
    if status != 0:
        print(f"set_default_handler: Launch Services error {status}", file=sys.stderr)
        return 1
    print(f"{content_type} files open with {bundle_id}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
