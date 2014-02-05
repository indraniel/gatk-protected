#!/bin/sh

mvn_args="verify"
mvn_properties=
mvn_clean=
unknown_args=
property_regex='-D(.*)=(.*)'
unit_test_regex='.*UnitTest'
post_script=
run_type="run"

for arg in "${@}" ; do
    if [[ "${arg}" == "dry" ]] ; then
        run_type="dry"

    elif [[ "${arg}" == "clean" ]] ; then
        mvn_clean="clean"
        mvn_args=

    elif [[ "${arg}" =~ ${property_regex} ]] ; then
        property_name=${BASH_REMATCH[1]}
        property_value=${BASH_REMATCH[2]}

        if [[ "${property_name}" == "single" ]] ; then
            test_property="test"
            test_disabled="it.test"
            if [[ ! "${property_value}" =~ ${unit_test_regex} ]] ; then
                test_property="it.test"
                test_disabled="test"
            fi

            mvn_properties="${mvn_properties} -D${test_disabled}=disabled -D${test_property}=${property_value}"

        elif [[ "${property_name}" == "test.debug.port" ]] ; then
            mvn_properties="${mvn_properties} -Dmaven.surefire.debug=\"-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=y,address=${property_value}\""
            mvn_properties="${mvn_properties} -Dmaven.failsafe.debug=\"-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=y,address=${property_value}\""

        elif [[ "${property_name}" == "test.default.maxmemory" ]] ; then
            mvn_properties="${mvn_properties} -Dtest.maxmemory=${property_value}"

        else
            unknown_args="${unknown_args} \"${arg}\""

        fi

    else
        if [[ "${arg}" != "dist" && "${mvn_args}" != "" && "${mvn_args}" != "verify" ]] ; then
            echo "Sorry, this script does not currently support mixing targets." >&2
            exit 1

        elif [[ "${arg}" == "dist" ]] ; then
            mvn_args="verify"

        elif [[ "${arg}" == "gatk" ]] ; then
            mvn_args="verify '-P!queue'"

        elif [[ "${arg}" == "test.compile" ]] ; then
            mvn_args="test-compile"

        elif [[ "${arg}" == "gatkdocs" ]] ; then
            local_repo="sitetemprepo"
            mvn_args="install -Dmaven.repo.local=${local_repo} -Ddisable.queue && mvn site -Dmaven.repo.local=${local_repo} -Ddisable.queue"

        elif [[ "${arg}" == "package.gatk.full" ]] ; then
            mvn_args="package '-P!private,!queue'"

        elif [[ "${arg}" == "package.gatk.all" ]] ; then
            mvn_args="package '-P!queue'"

        elif [[ "${arg}" == "package.queue.full" ]] ; then
            mvn_args="package '-P!private'"

        elif [[ "${arg}" == "package.queue.all" ]] ; then
            mvn_args="package"

#        elif [[ "${arg}" == "release.gatk.full" ]] ; then
#            mvn_args="package '-P!private,!queue'"
#            post_script=" && private/src/main/scripts/shell/copy_release.sh public/gatk-package/target/GenomeAnalysisTK-*.tar.bz2"

#        elif [[ "${arg}" == "release.queue.full" ]] ; then
#            mvn_args="package '-P!private'"
#            post_script=" && private/src/main/scripts/shell/copy_release.sh public/queue-package/target/Queue-*.tar.bz2"

        elif [[ "${arg}" == "build-picard-private" ]] ; then
            mvn_args="mvn install -f private/picard-maven/pom.xml"

        # TODO: clover support
        # see ant and maven docs for clover:
        #  https://confluence.atlassian.com/display/CLOVER/1.+QuickStart+Guide
        #  https://confluence.atlassian.com/display/CLOVER/Clover-for-Maven+2+and+3+User%27s+Guide
        #
        #elif [[ "${arg}" == "clover.report" ]] ; then
        #    mvn_args=...
        #
        #elif [[ "${arg}" == "with.clover" ]] ; then
        #    mvn_args=...

        # TODO: This runs *all* commit tests, including the few on Queue.
        elif [[ "${arg}" == "gatkfull.binary.release.tests" ]] ; then
            local_repo="sitetemprepo"
            mvn_args="install -Dmaven.repo.local=${local_repo} && mvn verify"
            mvn_args="${mvn_args} -Dmaven.repo.local=${local_repo}"
            mvn_args="${mvn_args} -Dsting.packagetests.enabled=true"
            mvn_args="${mvn_args} -Dsting.packagecommittests.skipped=false"

        # TODO: This runs only the pipeline tests (full, non-dry run), but not the commit tests for Queue.
        elif [[ "${arg}" == "queuefull.binary.release.tests" ]] ; then
            local_repo="sitetemprepo"
            mvn_args="install -Dmaven.repo.local=${local_repo} && mvn verify"
            mvn_args="${mvn_args} -Dmaven.repo.local=${local_repo}"
            mvn_args="${mvn_args} -Dsting.packagetests.enabled=true"
            mvn_args="${mvn_args} -Dsting.packagepipelinetests.skipped=false"
            mvn_args="${mvn_args} -Dsting.pipelinetests.run=true"

        elif [[ "${arg}" == "committests" ]] ; then
            mvn_args="verify -Dsting.committests.skipped=false"

        elif [[ "${arg}" == "test" ]] ; then
            mvn_args="test -Dsting.unittests.skipped=false"

        elif [[ "${arg}" == "unittest" ]] ; then
            mvn_args="test -Dsting.unittests.skipped=false"

        elif [[ "${arg}" == "integrationtest" ]] ; then
            mvn_args="verify -Dsting.integrationtests.skipped=false"

        elif [[ "${arg}" == "largescaletest" ]] ; then
            mvn_args="verify -Dsting.largescaletests.skipped=false"

        elif [[ "${arg}" == "pipelinetest" ]] ; then
            mvn_args="verify -Dsting.pipelinetests.skipped=false"

        elif [[ "${arg}" == "pipelinetestrun" ]] ; then
            mvn_args="verify -Dsting.pipelinetests.skipped=false -Dsting.pipelinetests.run=true"

        elif [[ "${arg}" == "fasttest" ]] ; then
            mvn_args="verify -Dsting.committests.skipped=false -pl private/gatk-private -am -Dresource.bundle.skip=true"

        else
            unknown_args="${unknown_args} \"${arg}\""

        fi

    fi

done

mvn_cmd=
if [[ "${mvn_clean}" != "" ]] ; then
    if [[ "${mvn_args}" != "" ]] ; then
      mvn_cmd="mvn ${mvn_clean} && mvn ${mvn_args}"
    else
      mvn_cmd="mvn ${mvn_clean}"
    fi
else
    mvn_cmd="mvn ${mvn_args}"
fi

if [[ "${unknown_args}" != "" ]] ; then
    echo "Unrecognized arguments:${unknown_args}" >&2

else
    echo "Equivalent maven command"
    echo "${mvn_cmd}${mvn_properties}${post_script}"

    if [[ "${run_type}" != "dry" ]] ; then
        sh -c "${mvn_cmd}${mvn_properties}${post_script}"
    fi

fi
