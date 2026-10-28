# AWS Batch dev environment

`aws-batch-dev.cfn.yml` sets up a small, dev-only AWS Batch environment for
running the test data and a few small samples. It isn't meant for production.

| Resource | Name |
| --- | --- |
| S3 bucket (work dirs, test data, results; `work/` expires after 14 days) | `atx-fastq2frags-dev-<account id>` |
| Batch compute environment: EC2 on-demand, m6i/r6i/c6i, 0–64 vCPU, 200 GB disk | `atx-fastq2frags-dev` |
| Batch job queue | `atx-fastq2frags-dev` |
| Instance role (ECS, SSM Session Manager, bucket read/write) | `atx-fastq2frags-batch-instance` |
| Seqera role, created only when `SeqeraExternalId` is set | `atx-fastq2frags-seqera` |

The compute environment scales to zero when nothing is queued, so it costs
nothing at idle. You pay for the bucket and for instances only while jobs run.

The commands below assume a shell with admin credentials for the target
account (aws-nonprod, 877829285799) in us-east-1.

## 1. Publish the container

Follow [`.github/aws/README.md`](../.github/aws/README.md) to deploy the ECR
stack and push a tag. The workflow publishes two images:

- `<ecr>/atx-fastq2frags:<version>`: the release build, with a 10 GiB Chromap sort buffer.
- `<ecr>/atx-fastq2frags:<version>-local`: a 1 GiB buffer. Use it with `-profile test`.

## 2. Deploy the Batch stack

```bash
aws cloudformation deploy \
  --stack-name fastq2frags-batch-dev \
  --template-file infra/aws-batch-dev.cfn.yml \
  --capabilities CAPABILITY_NAMED_IAM \
  --parameter-overrides \
      VpcId=vpc-02bce72d54b5c7bdf \
      SubnetIds=subnet-09122eeb006650918,subnet-03d2d742a6c58cde2,subnet-04858255ddcdf44c6
```

These are the public subnets of aws-nonprod's Control Tower VPC. The account
has no NAT gateway, so instances reach ECR and the internet through public IPs.
They have no inbound access. S3 traffic goes through the VPC's gateway endpoint.

## 3. Upload the test data

```bash
B=atx-fastq2frags-dev-877829285799
aws s3 sync tests/test_data/ s3://$B/test-data/ --exclude samplesheet.csv
cat > /tmp/samplesheet.s3.csv <<EOF
run_id,r1,r2,ng_id,ds_table
demo,s3://$B/test-data/demo_R1.fastq.gz,s3://$B/test-data/demo_R2.fastq.gz,,
demoxt,s3://$B/test-data/demo_R1.fastq.gz,s3://$B/test-data/demo_R2.fastq.gz,NG00001,s3://$B/test-data/xtalkTable.csv
EOF
aws s3 cp /tmp/samplesheet.s3.csv s3://$B/test-data/samplesheet.csv
```

## 4. Run from a laptop

Nextflow runs locally and submits each task to Batch. The laptop has to stay
online until the run finishes.

```bash
B=atx-fastq2frags-dev-877829285799
IMG=877829285799.dkr.ecr.us-east-1.amazonaws.com/atx-fastq2frags:$(cat version)-local

nextflow run main.nf -profile awsbatch,test \
    -w s3://$B/work \
    --container $IMG \
    --input s3://$B/test-data/samplesheet.csv --xtalk_correction \
    --genome_dir s3://$B/test-data/ref --barcode_file bc50.txt \
    --outdir s3://$B/results/test
```

The first task waits a few minutes while Batch starts an instance. Follow the
progress in the Nextflow output, in the Batch console, or in CloudWatch Logs
under `/aws/batch/job`.

## 5. Seqera Platform (trial)

Seqera uses role-based credentials, so no access keys are stored. Its Batch
Forge creates its own `TowerForge-*` queues and compute environments, and
deletes them when you remove the compute environment in Seqera.

1. In Seqera, go to **Credentials → Add → Amazon** and choose **role-based**.
   In **Assume role**, enter `arn:aws:iam::877829285799:role/atx-fastq2frags-seqera`.
   Seqera shows an **External ID**. Copy it, and save the credential after step 2.
2. Create or update the role with the External ID. Each workspace has its own
   ID; to use the role from several workspaces, pass them comma separated
   (`SeqeraExternalId=<id 1>,<id 2>`):

   ```bash
   aws cloudformation deploy --stack-name fastq2frags-batch-dev \
     --template-file infra/aws-batch-dev.cfn.yml --capabilities CAPABILITY_NAMED_IAM \
     --parameter-overrides SeqeraExternalId=<external id> \
         VpcId=vpc-02bce72d54b5c7bdf \
         SubnetIds=subnet-09122eeb006650918,subnet-03d2d742a6c58cde2,subnet-04858255ddcdf44c6
   ```

3. Go to **Compute environments → Add → AWS Batch → Batch Forge**:
   - Region: `us-east-1`
   - Work directory: `s3://atx-fastq2frags-dev-877829285799`
   - Provisioning: on-demand
   - Max CPUs: 128
   - Instance types: `m6i, r6i, c6i`
   - Boot disk size: 2000 GB. A real sample's FILTER_LINKERS needs about 2 TB
     of scratch (Latch allocates 1.95 TB).
   - Under advanced options, set the VPC and the three subnets above.

   Fusion v2 with fast instance storage also works (it needs Wave and the
   `atx-fastq2frags-wave-ecr` registry credential), but in our test it was
   slower and cost more. With instance types left empty, Batch picked older
   r5ad instances.
4. Add a GitHub credential, because the repo is private.
5. Add the pipeline from `https://github.com/atlasxomics/fastq2frags_nf` to the
   Launchpad. `nextflow_schema.json` provides the launch form.

## 6. Customer deliveries through Seqera (trial)

Customers get their results through Seqera's Data Explorer. Runs happen in an
internal workspace and publish to a delivery bucket per customer. Each customer
has a private workspace whose only credential is a role for their own bucket.

`delivery-customer.cfn.yml` (one stack per customer, `atx-delivery-<id>`) creates:

- `atx-delivery-<id>-877829285799`: SSE-S3, lifecycle expiry, CORS for Data
  Explorer, and a bucket policy that lets `TowerForge-*` instance roles publish
  and the internal workspace's credential (`atx-fastq2frags-seqera`) browse.
- `atx-delivery-<id>-seqera`: Seqera's documented minimum on that bucket only
  (list, read, upload), plus `s3:ListAllMyBuckets`, which Seqera needs to
  validate the credential.

A bucket per customer, rather than a prefix per customer in a shared bucket,
because Seqera's data-link check fails unless the role can list the whole
bucket, which would show customers each other's folder names.
`ListAllMyBuckets` still lets the role see other buckets' names (not contents);
Data Explorer only shows customers the buckets they can open.

For each customer:

1. In the customer's workspace, add role-based AWS credentials with
   `arn:aws:iam::877829285799:role/atx-delivery-<id>-seqera` and copy the External ID.
2. Deploy `delivery-customer.cfn.yml` with `CustomerId=<id>` and that External ID,
   then save the credential.
3. In **Data Explorer**, check that `s3://atx-delivery-<id>-877829285799` is
   listed (add it as a data link with that credential if not).
4. Invite the customer's users as collaborators with the **View** role.
5. Launch runs from the internal workspace with
   `--outdir s3://atx-delivery-<id>-877829285799/<run>`.

Keep that role as the only credential in the customer's workspace. Data
Explorer lists every bucket any workspace credential can read, and hiding a
data link doesn't block access.

Buckets are kept when a stack is deleted (`DeletionPolicy: Retain`); empty and
delete them in S3.

## Tear down

```bash
aws cloudformation delete-stack --stack-name fastq2frags-batch-dev
```

The bucket is retained when the stack is deleted. Empty it and delete it
separately. In Seqera, delete the compute environment to remove the
`TowerForge-*` resources.
