#!/usr/bin/env bash

set -euo pipefail
# --- CONFIGURATION ---
PROFILE="cnh-sso"

usage() {
  echo "Usage: $0 [BUCKET]"
  echo ""
  echo "Arguments:"
  echo ""
  echo "Example:"
  echo "  $0 bti-private-us-east-1-prd-gilbert-lab"
  echo ""
  exit 1
}

# Help argument
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
fi

# Parse argument for BUCKET
BUCKET="${1:-}"
if [[ -z "$BUCKET" ]]; then
  usage
fi

MOUNT_PATH="./data/$BUCKET"

echo "Using BUCKET: $BUCKET"


# --- AWS SSO Credentials for mount S3 ---

if ! aws sts get-caller-identity --profile "$PROFILE" > /dev/null 2>&1; then
  echo "🔐 AWS SSO session expired or missing. Logging in..."
  aws sso login --profile "$PROFILE"
  aws sts get-caller-identity --profile "$PROFILE" > /dev/null
else
  echo "✅ AWS SSO session is active."
fi

# check credentials file
CRED_FILE=$(ls -t ~/.aws/cli/cache/*.json | head -n 1)
if [[ ! -f "$CRED_FILE" ]]; then
  echo "❌ No valid AWS SSO credentials found. Please run 'aws sso login --profile $PROFILE' first."
else
  echo "✅ Found AWS SSO credentials file: $CRED_FILE"
fi

# export credentials
export AWS_ACCESS_KEY_ID=$(jq -r '.Credentials.AccessKeyId' < "$CRED_FILE")
export AWS_SECRET_ACCESS_KEY=$(jq -r '.Credentials.SecretAccessKey' < "$CRED_FILE")
export AWS_SESSION_TOKEN=$(jq -r '.Credentials.SessionToken' < "$CRED_FILE")

# --- Mount S3 Bucket ---
# Forcefully remove the mount point if it exists
if mountpoint -q "$MOUNT_PATH"; then
  echo "Forcefully removing existing mount point at $MOUNT_PATH..."
  sudo fuser -km "$MOUNT_PATH" 2>/dev/null || true
  fusermount -u "$MOUNT_PATH" 2>/dev/null || sudo umount -l "$MOUNT_PATH" 2>/dev/null || true
  rmdir "$MOUNT_PATH" 2>/dev/null || true
fi

# Create the mount directory
mkdir -p "$MOUNT_PATH"

# Mount the S3 bucket
echo "🔗 Mounting s3://$BUCKET to $MOUNT_PATH ..."
mount-s3 --allow-other "$BUCKET" "$MOUNT_PATH"

echo "🎉 Data and references preparation completed successfully."
echo "🎉 S3 bucket mount completed successfully."