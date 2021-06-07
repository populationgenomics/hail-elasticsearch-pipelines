version development

workflow hello {
  input {
    String? inp = "Hello, world!"
  }
  call echo as hello {
    input:
      inp=select_first([inp, "Hello, world!"])
  }
  output {
    File out = hello.out
  }
}

task echo {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    String inp
    Boolean? include_newline
  }
  command <<<
    set -e
    echo \
      ~{if (defined(include_newline) && select_first([include_newline])) then "-n" else ""} \
      '~{inp}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
    memory: "~{select_first([runtime_memory, 4])}G"
  }
  output {
    File out = stdout()
  }
}