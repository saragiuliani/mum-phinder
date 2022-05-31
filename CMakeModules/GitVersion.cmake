
if(__GitVersion)
    return()
endif()
set(__GitVersion YES)

get_filename_component(_GitVersionDir ${CMAKE_CURRENT_LIST_FILE} PATH)


function(SetGitVersionVariables)

    execute_process(COMMAND git log --pretty=format:'%h' -n 1
                    OUTPUT_VARIABLE GIT_REV
                    ERROR_QUIET)

    # Check whether we got any revision (which isn't
    # always the case, e.g. when someone downloaded a zip
    # file from Github instead of a checkout)
    if ("${GIT_REV}" STREQUAL "")
        set(GIT_VERSION_REV "N/A" PARENT_SCOPE)
        set(GIT_VERSION_DIFF "" PARENT_SCOPE)
        set(GIT_VERSION_TAG "N/A" PARENT_SCOPE)
        set(GIT_VERSION_BRANCH "N/A" PARENT_SCOPE)
        set(GIT_VERSION_DESC "N/A" PARENT_SCOPE)
    else()
        execute_process(
            COMMAND bash -c "git diff --quiet --exit-code || echo +"
            OUTPUT_VARIABLE GIT_DIFF)
        execute_process(
            COMMAND git rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE GIT_BRANCH)
        execute_process(
            COMMAND git rev-parse HEAD
            OUTPUT_VARIABLE GIT_COMMIT)
        execute_process(
            COMMAND git describe --tags
            OUTPUT_VARIABLE GIT_DESC
            ERROR_QUIET)

        string(STRIP "${GIT_REV}" GIT_REV)
        string(SUBSTRING "${GIT_REV}" 1 7 GIT_REV)
        string(STRIP "${GIT_DIFF}" GIT_DIFF)
        string(STRIP "${GIT_TAG}" GIT_TAG)
        string(STRIP "${GIT_BRANCH}" GIT_BRANCH)

        if ("${GIT_TAG}" STREQUAL "")
            set(GIT_TAG "v0.0.0")
        endif()

        if ("${GIT_DESC}" STREQUAL "")
            set(GIT_DESC "${GIT_TAG}-${GIT_REV}${GIT_DIFF}")
        endif()

        set(GIT_VERSION_REV "${GIT_REV}" PARENT_SCOPE)
        set(GIT_VERSION_DIFF "${GIT_DIFF}" PARENT_SCOPE)
        set(GIT_VERSION_TAG "${GIT_TAG}" PARENT_SCOPE)
        set(GIT_VERSION_DESC "${GIT_DESC}" PARENT_SCOPE)
        set(GIT_VERSION_BRANCH "${GIT_BRANCH}" PARENT_SCOPE)
    endif()

endfunction()