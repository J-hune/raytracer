execute_process(
    COMMAND "${BIN2C}" -c --name deviceCode "${INPUT}"
    OUTPUT_FILE "${OUTPUT}"
    COMMAND_ERROR_IS_FATAL ANY
)
file(READ "${OUTPUT}" source)
string(REPLACE "const unsigned char deviceCode[]" "extern const unsigned char deviceCode[]" source "${source}")
file(WRITE "${OUTPUT}" "${source}")
file(APPEND "${OUTPUT}" "\nextern \"C\" const unsigned long deviceCodeSize = sizeof(deviceCode);\n")
