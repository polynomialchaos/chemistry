################################################################################
# @file parse_chemkin.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
def read_encoded_file(path, encodings=None):
    """Return list of strings and try different encodings."""
    encodings = ['utf-8', 'latin-1'] if encodings is None else encodings

    for key in encodings:
        try:
            with open(path, 'r', encoding=key) as fptr:
                return fptr.readlines()[:]
        except UnicodeDecodeError:
            continue

    raise UnicodeDecodeError('Could not decode file: "{:}"!'.format(path))


def chemkin_format_reader(path, start_keys=None, end_keys=None, comment='!'):
    """Return list of strings and comments for a given path,
    start- and end-keys."""
    # sort keys by length to ensure the right replacement
    start_keys = sorted(
        [] if start_keys is None else start_keys, key=len, reverse=True)
    end_keys = sorted([] if end_keys is None else end_keys,
                      key=len, reverse=True)

    # read the file and check for encodings
    strings = read_encoded_file(path)

    # depending on the given start keys set the append mode
    app_mode = True if not start_keys else False

    # cycle through all lines and parse strings and comments
    parsed_data = []
    for _, line in enumerate(strings):
        # cycle commented or empty lines
        if not line.strip():
            continue
        if line.strip()[0] == comment:
            continue

        # start the append mode if a startKey is present in the current line
        if any(key in line for key in start_keys):
            app_mode = True

        # split the line into string and comment (keep comments might be needed)
        if app_mode and line.rstrip():
            tmp_line = line.split(comment, 1)
            parsed_data.append((
                tmp_line[0].rstrip(),
                tmp_line[1].rstrip() if tmp_line[1:] else ''
            ))

        # end append mode if no endKey is in strings or
        # if a endKey is present in the current line
        if app_mode and any(key in line for key in end_keys):
            break

    if parsed_data:
        # look for the keywords in the first line
        for key in start_keys:
            if key in parsed_data[0][0]:
                parsed_data[0] = (parsed_data[0][0].replace(
                    key, ''), parsed_data[0][1])

        # look for the keywords in the last line
        for key in end_keys:
            if key in parsed_data[-1][0]:
                parsed_data[-1] = (parsed_data[-1]
                                   [0].replace(key, ''), parsed_data[0][1])

    # clean empty lines
    parsed_data = [x for x in parsed_data if x[0].strip()]

    return parsed_data
