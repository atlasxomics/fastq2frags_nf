# AWS setup for the "Publish container" workflow

One-time setup, done by someone with IAM and ECR admin rights in the target
account. `ecr-publish.cfn.yml` creates:

- the ECR repository `atx-fastq2frags` (immutable tags, scan on push);
- the GitHub OIDC identity provider, unless it already exists;
- the IAM role `atx-fastq2frags-github-publish`, which only this GitHub
  repository can assume and which can only push to that ECR repository.

No AWS keys are stored in GitHub.

## 1. Deploy the stack

Check whether the account already has the GitHub OIDC provider:

```bash
aws iam list-open-id-connect-providers
```

If an entry ends in `token.actions.githubusercontent.com`, add
`CreateOIDCProvider=false` to `--parameter-overrides` below.

```bash
aws cloudformation deploy \
  --region us-east-1 \
  --stack-name fastq2frags-ecr-publish \
  --template-file .github/aws/ecr-publish.cfn.yml \
  --capabilities CAPABILITY_NAMED_IAM
```

Get the role ARN:

```bash
aws cloudformation describe-stacks \
  --region us-east-1 \
  --stack-name fastq2frags-ecr-publish \
  --query "Stacks[0].Outputs[?OutputKey=='RoleArn'].OutputValue" \
  --output text
```

## 2. Configure the GitHub repository

```bash
gh variable set AWS_ROLE_ARN --repo atlasxomics/fastq2frags_nf --body "<role arn>"
```

Set `AWS_REGION` or `ECR_REPOSITORY` the same way only if you changed them
from `us-east-1` / `atx-fastq2frags`.

## 3. Publish

```bash
git tag "v$(cat version)"
git push origin "v$(cat version)"
```
