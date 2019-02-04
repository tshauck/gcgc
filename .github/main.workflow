workflow "New workflow" {
  on = "push"
  resolves = ["Build Image", "Run Tests"]
}

action "Build Image" {
  uses = "actions/docker/cli@c08a5fc9e0286844156fefff2c141072048141f6"
  args = "build -t thauck/gcgc-tester ."
}

action "Run Tests" {
  uses = "actions/docker/cli@c08a5fc9e0286844156fefff2c141072048141f6"
  args = "run --rm -w=/gcgc --entrypoint=pytest -t thauck/gcgc-tester"
  needs = ["Build Image"]
}
