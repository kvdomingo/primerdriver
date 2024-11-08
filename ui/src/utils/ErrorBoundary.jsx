import { MDBTypography as Typography } from "mdbreact";
import { Component } from "react";

class ErrorBoundary extends Component {
  state = {
    hasError: false,
  };

  static getDerivedStateFromError(_) {
    return { hasError: true };
  }

  componentDidCatch(error, errorInfo) {}

  render() {
    return this.state.hasError ? (
      <Typography
        variant="h5-responsive"
        tag="p"
        className="text-center align-items-center"
      >
        Something went wrong 😢
      </Typography>
    ) : (
      this.props.children
    );
  }
}

export default ErrorBoundary;
