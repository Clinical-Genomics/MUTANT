### The purpose of the code changes are as follows:
-  THING

### Preparations:

**How to prepare mutant for test**:
* `bash /home/proj/stage/mutant/MUTANT/mutant/standalone/deploy_hasta_update.sh stage <MUTANT_branch>`

**How to prepare cg for test**:
- `us`
- Paxa and update cg or other tools needed for the test:
  - `paxa -u <user> -s hasta -r cg-stage`
  - `paxa -u <user> -s hasta -r servers-stage`
  - `bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-cg-stage.sh master`
  - `bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-servers-stage.sh master`
  

### How to test:
Run in a tmux screen or similar if testing on a large dataset.

**Test with cg**:
- `us`
- `cg workflow mutant start maturejay`

### Expected outcome:
- [ ] Produced files contain expected values

### Review:
- [ ] Code reviewed by

This [version](https://semver.org/) is a:
- [ ] **MAJOR** - when you make incompatible API changes
- [ ] **MINOR** - when you add functionality in a backwards compatible manner
- [ ] **PATCH** - when you make backwards compatible bug fixes or documentation/instructions
