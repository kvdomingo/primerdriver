import { MDBTypography as Typography } from "mdbreact";
import { Component } from "react";

class ErrorBoundary extends Component {
  state = {
    hasError: false,
  };

  static getDerivedStateFromError(error) {
    return { hasError: true };
  }

  componentDidCatch(error, errorInfo) {}

  render() {
    return this.state.hasError ? (
      <Typography variant="h5-responsive" tag="p" className="align-items-center text-center">
        Something went wrong 😢
      </Typography>
    ) : (
      this.props.children
    );
  }
}

export default ErrorBoundary;
